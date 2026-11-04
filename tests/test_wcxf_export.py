"""Check operator coverage, conventions and units at the Wilson/WCxf boundary."""

import numpy as np
import pytest
from wilson import Wilson, wcxf

from lsme.numeric import GranadaNMatchingResult as N
from lsme.numeric import GranadaTheta1MatchingResult as Theta1
from lsme.numeric import GranadazetaMatchingResult as Zeta
from lsme.numeric.matchingresult import GenericMatchingResult


def test_export_contains_all_baryon_conserving_operators():
    sectors = wcxf.Basis["SMEFT", "Warsaw"].sectors
    baryon_violating_sectors = {"dB=de=1", "dB=dmu=1", "dB=dtau=1", "dB=dL=1"}
    expected = {
        op
        for sector, operators in sectors.items()
        if sector not in baryon_violating_sectors
        for op in operators
    }
    result = Zeta()
    coefficients = result.coeff_dict()
    wc = result.to_wcxf()

    assert set(coefficients) == expected
    assert set(wc.dict) == expected
    assert "llphiphi_11" in coefficients
    assert (wc.eft, wc.basis) == ("SMEFT", "Warsaw")
    wc.validate()


def test_export_does_not_depend_on_sector_names(monkeypatch):
    basis = wcxf.Basis["SMEFT", "Warsaw"]
    result = Zeta()
    expected = result.coeff_dict()
    expected_wcxf = result.to_wcxf().dict
    # Put supported and unsupported operators together under an unfamiliar
    # label: exclusion must follow the operators rather than their grouping.
    combined = {op: metadata for operators in basis.sectors.values() for op, metadata in operators.items()}
    monkeypatch.setattr(basis, "sectors", {"regrouped": combined})

    assert result.coeff_dict() == expected
    assert result.to_wcxf().dict == expected_wcxf


def test_dimension_six_coefficients_use_tev_and_wcxf_uses_gev():
    result = Theta1()
    result.onelooporder = 0
    result.MTheta1 = 2.0
    result.scale = 1.5
    result.lambdaTheta1 = 0.6 + 0.2j
    result.lambdaTheta1bar = result.lambdaTheta1.conjugate()

    expected_tev = abs(result.lambdaTheta1)**2 / (6 * 2.0**2)
    expected_gev = abs(result.lambdaTheta1)**2 / (6 * 2000.0**2)
    assert result.alphaOH() == pytest.approx(expected_tev, rel=1e-13, abs=0)
    assert result.coeff_dict()["phi"] == result.alphaOH()

    wc = result.to_wcxf()
    wc.validate()
    assert wc.scale == pytest.approx(1500.0, rel=1e-15)
    assert wc.dict["phi"] == pytest.approx(expected_gev, rel=1e-13, abs=0)
    assert result.scale == 1.5

    w = Wilson.from_wc(wc)
    assert w.wc["phi"] == pytest.approx(expected_gev, rel=1e-13, abs=0)


def test_weinberg_coefficients_use_tev_and_wcxf_uses_gev_with_flavor_and_phase():
    result = N()
    result.onelooporder = 0
    result.MN = 2.0
    result.lambdaN = np.array([0.2 + 0.1j, 0.3 - 0.4j, 0.5 + 0.2j])
    result.lambdaNbar = result.lambdaN.conjugate()
    coefficients = result.coeff_dict()
    wc = result.to_wcxf()
    wc.validate()
    w = Wilson.from_wc(wc)

    for i in range(3):
        for j in range(i, 3):
            op = f"llphiphi_{i + 1}{j + 1}"
            expected_tev = result.lambdaN[i] * result.lambdaN[j] / (2 * 2.0)
            assert result.alphaWeinberg(i, j) == pytest.approx(expected_tev, rel=1e-13, abs=0)
            assert coefficients[op] == pytest.approx(expected_tev, rel=1e-13, abs=0)
            assert wc.dict[op] == pytest.approx(expected_tev / 1000.0, rel=1e-13, abs=0)
            assert w.wc[op] == pytest.approx(expected_tev / 1000.0, rel=1e-13, abs=0)


def test_missing_supported_coefficient_is_not_silently_zero():
    incomplete = GenericMatchingResult("Incomplete")
    with pytest.raises(AttributeError, match="alphaO3G"):
        incomplete.coeff_dict()
    with pytest.raises(AttributeError, match="alphaO3G"):
        incomplete.to_wcxf()


@pytest.mark.parametrize("operator", ["ll", "qq1", "qq3", "uu", "dd", "ee"])
def test_full_current_tensor_is_collected_before_export(operator, monkeypatch):
    # Hermitian coefficients with deliberately unequal equivalent entries.
    # An orbit sum is an independent Lagrangian-level oracle; multiplying a
    # single entry by its multiplicity fails this test.
    rng = np.random.default_rng(12)
    raw = rng.normal(size=(3,) * 4) + 1j * rng.normal(size=(3,) * 4)
    raw = (raw + raw.transpose(1, 0, 3, 2).conjugate()) / 2
    result = Zeta()
    method = result._warsaw_method(operator)[0]
    monkeypatch.setattr(result, method, lambda *indices: raw[indices])
    native = result.coeff_dict()
    exported = result.to_wcxf().dict
    multiplicities = set()
    for key in native:
        if not key.startswith(operator + "_"):
            continue
        p, r, s, t = (int(i) - 1 for i in key.split("_")[1])
        orbit = {(p, r, s, t), (s, t, p, r)}
        if operator == "ee":
            orbit |= {(s, r, p, t), (p, t, s, r)}
        multiplicities.add(len(orbit))
        assert native[key] == raw[p, r, s, t]
        assert exported[key] == pytest.approx(sum(raw[i] for i in orbit) * 1e-6)
    assert multiplicities == ({1, 2, 4} if operator == "ee" else {1, 2})
    assert result.coeff_dict() == native


def test_export_projects_only_real_warsaw_keys(monkeypatch):
    result = Theta1()
    raw = 2 + 3j
    for family in ("phie", "ee", "llphiphi"):
        method, _ = result._warsaw_method(family)
        monkeypatch.setattr(result, method, lambda *indices: raw)

    values = result.to_wcxf().dict
    assert values["phie_11"] == pytest.approx(raw.real * 1e-6)
    assert values["phie_12"] == pytest.approx(raw * 1e-6)
    assert values["ee_1122"] == pytest.approx(4 * raw.real * 1e-6)
    assert values["llphiphi_11"] == pytest.approx(raw * 1e-3)


def test_s2_tree_ee_has_factor_four_but_e_gauge_loop_has_factor_two():
    from lsme.numeric import GranadaS2MatchingResult as S2
    from lsme.numeric import GranadaEMatchingResult as E

    scalar = S2()
    scalar.onelooporder = 0
    assert scalar.coeff_dict()["ee_1122"] == pytest.approx(0.5)
    assert scalar.to_wcxf().dict["ee_1122"] == pytest.approx(2e-6)

    # E's pure-gauge ee tensor is a delta_pr delta_st, not Fierz-symmetric.
    # Only 1122 and 2211 contribute, even though the orbit has four entries.
    lepton = E()
    a = -lepton.g1**4 / (120 * np.pi**2 * lepton.ME**2)
    assert lepton.alphaOee(0, 1, 1, 0) == 0
    assert lepton.coeff_dict()["ee_1122"] == pytest.approx(a)
    assert lepton.to_wcxf().dict["ee_1122"] == pytest.approx(2 * a * 1e-6)


@pytest.mark.parametrize("operator", [
    "uB", "uW", "uG", "dB", "dW", "dG", "eB", "eW",
    "G", "Gtilde", "W", "Wtilde",
    "phiB", "phiW", "phiG", "phiWB", "phiBtilde", "phiWtilde",
    "phiGtilde", "phiWtildeB",
])
def test_field_strength_signs_and_complex_phases(operator, monkeypatch):
    result = Zeta()
    dipole = operator[0] in "ude"
    raw = 0.7 + 0.2j if dipole else 0.7
    monkeypatch.setattr(result, result._warsaw_method(operator)[0], lambda *indices: raw)
    native = result.coeff_dict()
    wc = result.to_wcxf()
    wc.validate()
    sign = 1 if operator.startswith("phi") else -1
    keys = [k for k in native if k.split("_")[0] == operator]
    assert len(keys) == (9 if dipole else 1)
    for key in keys:
        assert native[key] == raw
        assert wc.dict[key] == pytest.approx(sign * raw * 1e-6)
