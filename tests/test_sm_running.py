"""Check the map from Wilson's SM parameters to LSME and its use by the classes."""

import time

import numpy as np
import pytest
import ckmutil.ckm
from wilson.run.smeft import smpar

from lsme.numeric import GranadaSMatchingResult as S
from lsme.numeric import GranadaXiMatchingResult as Xi
from lsme.numeric import matchingresult
from lsme.numeric.matchingresult import GEV, GenericMatchingResult
from lsme.numeric.sm_running import SM_PARAMETER_NAMES, sm_parameters, wilson_sm_lagrangian
from lsme.numeric.sm_table import MZ_GEV, MU_MAX_GEV, VALIDATION_RTOL

MZ_TEV = MZ_GEV * GEV
SCALES_TEV = [MZ_TEV, 0.25, 1.0, 7.3, 120.0, MU_MAX_GEV * GEV]


@pytest.mark.parametrize("scale", SCALES_TEV)
def test_map_from_wilson_is_exact(scale):
    lsme = sm_parameters(scale)
    wilson = wilson_sm_lagrangian(scale / GEV)

    assert set(lsme) == set(SM_PARAMETER_NAMES)
    assert (lsme["g1"], lsme["g2"], lsme["g3"]) == (wilson["gp"], wilson["g"], wilson["gs"])
    assert lsme["lam"] == wilson["Lambda"] / 2
    # MatchMaker: −muH^2 H†H; Wilson: +m2 φ†φ.
    assert lsme["muH"].real == 0
    assert lsme["muH"]**2 == pytest.approx(-wilson["m2"] * GEV**2, rel=1e-15, abs=0)
    for lsme_name, wilson_name in (("yu", "Gu"), ("yd", "Gd"), ("yl", "Ge")):
        np.testing.assert_array_equal(lsme[lsme_name], wilson[wilson_name])
        np.testing.assert_array_equal(lsme[lsme_name + "bar"], np.conjugate(wilson[wilson_name]))
        assert not np.shares_memory(lsme[lsme_name], lsme[lsme_name + "bar"])


def test_parameters_reproduce_wilson_mz_inputs():
    # The table approximates Wilson's extraction, so the reconstruction is
    # judged with the table's declared tolerance, not exactly.
    p = smpar.p
    lsme = sm_parameters(MZ_TEV)
    v = (np.sqrt(2) * p["GF"])**-0.5

    assert lsme["g3"] == pytest.approx(np.sqrt(4 * np.pi * p["alpha_s"]), rel=VALIDATION_RTOL)
    assert lsme["g2"] == pytest.approx(2 * p["m_W"] / v, rel=VALIDATION_RTOL)
    e = lsme["g1"] * lsme["g2"] / np.hypot(lsme["g1"], lsme["g2"])
    assert e**2 / (4 * np.pi) == pytest.approx(p["alpha_e"], rel=VALIDATION_RTOL)
    assert 2 * lsme["lam"] * v**2 == pytest.approx(p["m_h"]**2, rel=VALIDATION_RTOL)
    assert -lsme["muH"]**2 / GEV**2 == pytest.approx(p["m_h"]**2 / 2, rel=VALIDATION_RTOL)

    for name, masses in (
        ("yu", ("m_u", "m_c", "m_t")),
        ("yd", ("m_d", "m_s", "m_b")),
        ("yl", ("m_e", "m_mu", "m_tau")),
    ):
        singular = np.sort(np.linalg.svd(lsme[name], compute_uv=False))
        expected = np.array([p[m] for m in masses]) * np.sqrt(2) / v
        np.testing.assert_allclose(singular, expected, rtol=VALIDATION_RTOL)

    Uu, _, _ = np.linalg.svd(lsme["yu"])
    Ud, _, _ = np.linalg.svd(lsme["yd"])
    ckm = Uu[:, ::-1].conj().T @ Ud[:, ::-1]
    expected_ckm = ckmutil.ckm.ckm_tree(p["Vus"], p["Vub"], p["Vcb"], p["gamma"])
    np.testing.assert_allclose(np.abs(ckm), np.abs(expected_ckm), rtol=VALIDATION_RTOL)


def test_imaginary_muH_gives_real_coefficients():
    result = S()
    assert result.muH.real == 0 and result.muH.imag > 0
    assert (result.muH)**(2) == complex(-result.muH.imag**2, 0.0)

    # With real inputs the muH-dependent coefficients are real numbers held
    # as complex values with zero imaginary part.
    for name in ("yu", "yd", "yl"):
        real = np.diag(np.abs(np.diag(getattr(result, name))))
        setattr(result, name, real)
        setattr(result, name + "bar", real.copy())
    result.kappaS = result.kappaSbar = 0.4
    for method in ("alphaOH", "alphaOuH"):
        args = (2, 2) if method == "alphaOuH" else ()
        value = getattr(result, method)(*args)
        assert np.imag(value) == 0
        assert np.real(value) != 0


def test_exports_validate_with_imaginary_muH():
    result = Xi()
    assert isinstance(result.muH, complex)
    wc = result.to_wcxf()
    wc.validate()
    assert all(np.isfinite(complex(value)) for value in result.coeff_dict().values())


def test_scale_changes_couplings():
    result = GenericMatchingResult("Test")
    low = {name: getattr(result, name) for name in ("g3", "lam")}
    yt_low = np.abs(result.yu[2, 2])
    result.scale = 50.0
    assert result.scale == 50.0
    assert result.g3 < low["g3"] and result.lam < low["lam"]
    assert np.abs(result.yu[2, 2]) < yt_low

    again = GenericMatchingResult("Test", scale=50.0)
    for name in SM_PARAMETER_NAMES:
        np.testing.assert_array_equal(getattr(again, name), getattr(result, name))
        assert not np.shares_memory(getattr(again, name), getattr(result, name)) or np.ndim(getattr(result, name)) == 0


def _state(result):
    return {name: np.copy(getattr(result, name)) for name in ("scale",) + SM_PARAMETER_NAMES}


def _assert_same_state(result, state):
    for name, value in state.items():
        np.testing.assert_array_equal(getattr(result, name), value)


@pytest.mark.parametrize("scale", [np.nan, np.inf, -1.0, 0.0, 0.05, 1000.1])
def test_invalid_scale_leaves_object_unchanged(scale):
    result = GenericMatchingResult("Test", scale=2.0)
    state = _state(result)
    with pytest.raises(ValueError):
        result.scale = scale
    _assert_same_state(result, state)


def test_lookup_failure_leaves_object_unchanged(monkeypatch):
    result = GenericMatchingResult("Test", scale=2.0)
    state = _state(result)

    def failing_lookup(scale):
        raise RuntimeError("lookup failed")

    monkeypatch.setattr(matchingresult, "sm_parameters", failing_lookup)
    with pytest.raises(RuntimeError, match="lookup failed"):
        result.scale = 3.0
    _assert_same_state(result, state)


def test_just_running_freezes_couplings():
    result = S()
    result.MS = 3.0
    result.scale = 0.8
    result.kappaS = result.kappaSbar = 0.5
    state = _state(result)

    running = result.just_running("alphaOuH", [2, 2])
    _assert_same_state(result, state)

    at_scale = result.alphaOuH(2, 2)
    result._scale = result.MS
    at_mass = result.alphaOuH(2, 2)
    result._scale = 0.8
    assert running == at_scale - at_mass
    assert running != 0


def test_just_running_restores_scale_when_coefficient_raises():
    result = S()
    state = _state(result)
    result.alphaBroken = lambda: 1 / 0
    with pytest.raises(ZeroDivisionError):
        result.just_running("alphaBroken")
    _assert_same_state(result, state)


def test_just_running_accepts_explicit_reference_scale():
    result = S()
    result.kappaS = result.kappaSbar = 0.5
    expected = result.alphaOuH(2, 2)
    result._scale = 4.0
    expected -= result.alphaOuH(2, 2)
    result._scale = 1.0
    assert result.just_running("alphaOuH", [2, 2], reference_scale=4.0) == expected


@pytest.mark.parametrize("reference", [np.nan, np.inf, -np.inf, 0.0, -2.0])
def test_just_running_rejects_invalid_reference_scale(reference):
    result = S()
    state = _state(result)
    with pytest.raises(ValueError):
        result.just_running("alphaOuH", [2, 2], reference_scale=reference)
    _assert_same_state(result, state)


def test_just_running_accepts_reference_outside_table():
    # The reference scale only enters logarithms, so it need not lie in the
    # table's range.
    result = S()
    result.kappaS = result.kappaSbar = 0.5
    assert np.isfinite(result.just_running("alphaOuH", [2, 2], reference_scale=5000.0))


@pytest.mark.performance
def test_setting_scale_is_fast():
    result = GenericMatchingResult("Test")
    scales = np.geomspace(0.1, 900.0, 500)
    times = []
    for scale in scales:
        start = time.perf_counter()
        result.scale = scale
        times.append(time.perf_counter() - start)
    assert np.median(times) < 100e-6
