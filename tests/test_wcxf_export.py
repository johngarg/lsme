"""Check operator coverage and physical units at the Wilson/WCxf boundary."""

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
    # Put supported and unsupported operators together under an unfamiliar
    # label: exclusion must follow the operators rather than their grouping.
    combined = {op: metadata for operators in basis.sectors.values() for op, metadata in operators.items()}
    monkeypatch.setattr(basis, "sectors", {"regrouped": combined})

    assert result.coeff_dict() == expected


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
