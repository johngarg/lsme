import inspect
import json
from pathlib import Path

import numpy as np
import pytest

from lsme.numeric import GranadavarphiMatchingResult as Varphi


ARRAY_NAMES = (
    "yu",
    "yubar",
    "yd",
    "ydbar",
    "yl",
    "ylbar",
    "yvarphiu",
    "yvarphiubar",
    "yvarphid",
    "yvarphidbar",
    "yvarphie",
    "yvarphiebar",
)
BASELINE_PATH = Path(__file__).with_name("varphi_zero_coupling_baseline.json")


def configure_baseline_point(result, point):
    scalars = {
        "Mvarphi": 2.3 + 0.4 * point,
        "scale": 1.1 + 0.2 * point,
        "g1": 0.31 + 0.01 * point,
        "g2": 0.59 + 0.02 * point,
        "g3": 1.07 + 0.03 * point,
        "mH": 0.126 + 0.001 * point,
        "vev": 0.247 + 0.001 * point,
        "invepsilonbar": 0.13 + 0.02 * point,
        "epsilonbar": -0.07 + 0.01 * point,
        "onelooporder": 1,
        "lambdavarphi": 0.17 + 0.03 * point,
        "lambdavarphibar": -0.21 + 0.02 * point,
        "lambdaHatvarphi": 0.23 + 0.01 * point,
        "lambdaHatvarphibar": 0.19 + 0.01 * point,
        "lambdaHatPrimevarphi": -0.11 + 0.02 * point,
        "lambdaHatPrimevarphibar": 0.14 + 0.02 * point,
    }
    for name, value in scalars.items():
        setattr(result, name, value)
    for offset, name in enumerate(ARRAY_NAMES):
        value = (np.arange(9).reshape(3, 3) + 1 + offset + point) / (
            31.0 + offset
        )
        setattr(result, name, value)


def zero_exotic_yukawas(result):
    for name in (
        "yvarphiu",
        "yvarphiubar",
        "yvarphid",
        "yvarphidbar",
        "yvarphie",
        "yvarphiebar",
    ):
        setattr(result, name, np.zeros((3, 3)))


def coupling_shift(method, configure, args=()):
    with_coupling = Varphi()
    without_coupling = Varphi()
    configure(with_coupling)
    configure(without_coupling)
    without_coupling.lambdaHatPrimePrimevarphi = 0
    without_coupling.lambdaHatPrimePrimevarphibar = 0
    return (
        getattr(with_coupling, method)(*args)
        - getattr(without_coupling, method)(*args),
        with_coupling,
    )


def test_prime_prime_public_interface():
    result = Varphi()

    assert result.lambdaHatPrimePrimevarphi == 1
    assert result.lambdaHatPrimePrimevarphibar == 1
    assert "lambdaHatPrimePrimevarphi" in result.exotic_params
    assert "lambdaHatPrimePrimevarphibar" in result.exotic_params
    assert isinstance(result.coeff_dict(), dict)


@pytest.mark.parametrize("point", [0, 1])
def test_prime_prime_zero_reproduces_baseline(point):
    with BASELINE_PATH.open() as baseline_file:
        expected = json.load(baseline_file)[str(point)]

    result = Varphi()
    configure_baseline_point(result, point)
    result.lambdaHatPrimePrimevarphi = 0
    result.lambdaHatPrimePrimevarphibar = 0

    methods = sorted(
        name
        for name in dir(result)
        if name.startswith("alpha") and callable(getattr(result, name))
    )
    assert methods == sorted(expected)

    for name in methods:
        method = getattr(result, name)
        zero_flavour_indices = [0] * len(inspect.signature(method).parameters)
        assert method(*zero_flavour_indices) == pytest.approx(
            expected[name], rel=2e-12, abs=2e-12
        )


def configure_bosonic(result):
    result.Mvarphi = 2.7
    result.scale = 1.4
    result.lambdaHatPrimePrimevarphi = 0.37
    result.lambdaHatPrimePrimevarphibar = -0.29
    result.lambdaHatvarphi = 0.23
    result.lambdaHatPrimevarphi = -0.17
    result.lambdavarphi = 0
    result.lambdavarphibar = 0
    zero_exotic_yukawas(result)


def test_prime_prime_bosonic_shifts():
    x = 0.37
    xbar = -0.29

    box_shift, result = coupling_shift("alphaOHBox", configure_bosonic)
    hd_shift, _ = coupling_shift("alphaOHD", configure_bosonic)
    h_shift, _ = coupling_shift("alphaOH", configure_bosonic)

    common = x * xbar * result.onelooporder / (np.pi**2 * result.Mvarphi**2)
    assert box_shift == pytest.approx(common / 48)
    assert hd_shift == pytest.approx(common / 24)
    assert h_shift == pytest.approx(
        common
        * (
            2 * result.lam
            + 3 * result.lambdaHatvarphi
            + 3 * result.lambdaHatPrimevarphi
        )
        / 24
    )


def configure_phase_sensitive(result):
    result.Mvarphi = 3.1
    result.scale = 1.2
    result.invepsilonbar = 0.19
    result.lambdaHatPrimePrimevarphi = 0.41
    result.lambdaHatPrimePrimevarphibar = -0.27
    result.lambdavarphi = 0.13
    result.lambdavarphibar = -0.22
    result.lambdaHatvarphi = 0
    result.lambdaHatPrimevarphi = 0
    zero_exotic_yukawas(result)


def test_prime_prime_phase_sensitive_higgs_shift():
    shift, result = coupling_shift("alphaOH", configure_phase_sensitive)
    x = result.lambdaHatPrimePrimevarphi
    xbar = result.lambdaHatPrimePrimevarphibar
    product_term = (
        x
        * xbar
        * 2
        * result.lam
        * result.onelooporder
        / (24 * np.pi**2 * result.Mvarphi**2)
    )
    bracket = (
        5 * result.invepsilonbar
        + 7
        - 5 * np.log(result.Mvarphi**2 / result.mu**2)
    )
    expected = (
        (x * result.lambdavarphi**2 + xbar * result.lambdavarphibar**2)
        * bracket
        * result.onelooporder
        / (8 * np.pi**2 * result.Mvarphi**2)
    )

    assert shift - product_term == pytest.approx(expected)


@pytest.mark.parametrize(
    ("method", "sm_yukawa"),
    [("alphaOdH", "yd"), ("alphaOeH", "yl"), ("alphaOuH", "yu")],
)
def test_prime_prime_hbox_eom_propagation(method, sm_yukawa):
    def configure(result):
        configure_bosonic(result)
        setattr(
            result,
            sm_yukawa,
            (np.arange(9).reshape(3, 3) + 1) / 17.0,
        )

    shift, result = coupling_shift(method, configure, args=(1, 2))
    expected = (
        result.lambdaHatPrimePrimevarphi
        * result.lambdaHatPrimePrimevarphibar
        * getattr(result, sm_yukawa)[1, 2]
        * result.onelooporder
        / (48 * np.pi**2 * result.Mvarphi**2)
    )

    assert shift == pytest.approx(expected)


def configure_yukawa_dependent(result):
    result.Mvarphi = 2.9
    result.scale = 1.3
    result.invepsilonbar = 0.16
    result.lambdaHatPrimePrimevarphi = 0.31
    result.lambdaHatPrimePrimevarphibar = -0.28
    result.lambdavarphi = 0.21
    result.lambdavarphibar = -0.17
    for offset, name in enumerate(ARRAY_NAMES):
        setattr(
            result,
            name,
            (np.arange(9).reshape(3, 3) + 2 + offset) / (23.0 + offset),
        )


@pytest.mark.parametrize(
    ("method", "sm_name", "exotic_name", "coupling_name", "linear_name", "sign", "orientation"),
    [
        ("alphaOdH", "ydbar", "yvarphidbar", "lambdaHatPrimePrimevarphibar", "lambdavarphibar", 1, "down"),
        ("alphaOdHbar", "yd", "yvarphid", "lambdaHatPrimePrimevarphi", "lambdavarphi", 1, "down"),
        ("alphaOeH", "ylbar", "yvarphiebar", "lambdaHatPrimePrimevarphibar", "lambdavarphibar", 1, "down"),
        ("alphaOeHbar", "yl", "yvarphie", "lambdaHatPrimePrimevarphi", "lambdavarphi", 1, "down"),
        ("alphaOuH", "yubar", "yvarphiu", "lambdaHatPrimePrimevarphi", "lambdavarphi", -1, "up"),
        ("alphaOuHbar", "yu", "yvarphiubar", "lambdaHatPrimePrimevarphibar", "lambdavarphibar", -1, "up"),
    ],
)
def test_prime_prime_additional_exotic_yukawa_terms(
    method,
    sm_name,
    exotic_name,
    coupling_name,
    linear_name,
    sign,
    orientation,
):
    i, j = 1, 2
    shift, result = coupling_shift(method, configure_yukawa_dependent, args=(i, j))
    x = result.lambdaHatPrimePrimevarphi
    xbar = result.lambdaHatPrimePrimevarphibar
    sm_eom_name = {
        "alphaOdH": "yd",
        "alphaOdHbar": "ydbar",
        "alphaOeH": "yl",
        "alphaOeHbar": "ylbar",
        "alphaOuH": "yu",
        "alphaOuHbar": "yubar",
    }[method]
    eom_term = (
        x
        * xbar
        * getattr(result, sm_eom_name)[i, j]
        * result.onelooporder
        / (48 * np.pi**2 * result.Mvarphi**2)
    )

    sm = getattr(result, sm_name)
    exotic = getattr(result, exotic_name)
    if orientation == "down":
        cubic_yukawa = sum(
            sm[a, b] * exotic[b, i] * exotic[j, a]
            for a in range(3)
            for b in range(3)
        )
        linear_yukawa = exotic[j, i]
    else:
        cubic_yukawa = sum(
            sm[b, a] * exotic[b, j] * exotic[i, a]
            for a in range(3)
            for b in range(3)
        )
        linear_yukawa = exotic[i, j]

    bracket = 1 + result.invepsilonbar - np.log(result.Mvarphi**2 / result.mu**2)
    expected_extra = (
        getattr(result, coupling_name)
        * (
            cubic_yukawa
            + sign * 5 * getattr(result, linear_name) * bracket * linear_yukawa
        )
        * result.onelooporder
        / (8 * np.pi**2 * result.Mvarphi**2)
    )

    assert shift - eom_term == pytest.approx(expected_extra)
