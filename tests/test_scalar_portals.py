"""Scalar inputs use the potential convention of Table 2 of arXiv:2412.01759.

Tree references are from the Granada dictionary, arXiv:1711.10391, App. A.2
and D. For Xi1 we use the explicit halves on the two portals in Eq. (A.7).
"""

import numpy as np
import pytest

from lsme.numeric import (
    GranadaPhiMatchingResult,
    GranadaSMatchingResult,
    GranadaXiMatchingResult,
    GranadaXi1MatchingResult,
)


@pytest.mark.parametrize("mass", [1.0, 2.3])
def test_singlet_tree(mass):
    model = GranadaSMatchingResult()
    model.onelooporder = 0
    model.MS, model.kappaS, model.kappaS3, model.lambdaS = mass, 0.7, 0.2, 0.4

    # Eq. (D.42), including the cubic exchange with its relative sign.
    expected = 0.7**3 * 0.2 / mass**6 - 0.7**2 * 0.4 / mass**4
    assert model.alphaOH() == pytest.approx(expected)
    assert model.alphaOHBox() == pytest.approx(-0.7**2 / (2 * mass**4))
    assert model.coeff_dict()["phi"] == pytest.approx(expected)
    assert model.to_wcxf().dict["phi"] == pytest.approx(1e-6 * expected, abs=1e-18)


def test_real_triplet_tree():
    model = GranadaXiMatchingResult()
    model.onelooporder = 0
    model.MXi, model.kappaXi, model.lambdaXi = 1.8, -0.6, 0
    baseline = model.alphaOH()
    model.lambdaXi = -0.35
    k2_over_m4 = model.kappaXi**2 / model.MXi**4

    # Subtract the portal-independent tree terms, including the SM quartic.
    assert model.alphaOH() - baseline == pytest.approx(-model.lambdaXi * k2_over_m4)
    assert model.alphaOHD() == pytest.approx(-2 * k2_over_m4)


@pytest.mark.parametrize(("portal", "primed"), [(0.4, 0), (0, -0.3), (0.4, -0.3)])
def test_complex_triplet_tree(portal, primed):
    model = GranadaXi1MatchingResult()
    model.onelooporder = 0
    model.MXi1, model.kappaXi1 = 2.1, 0.4 + 0.3j
    model.kappaXi1bar = model.kappaXi1.conjugate()
    model.lambdaXi1 = model.lambdaXi1bar = 0
    model.lambdaXi1P = model.lambdaXi1Pbar = 0
    baseline = model.alphaOH()
    model.lambdaXi1 = model.lambdaXi1bar = portal
    model.lambdaXi1P = model.lambdaXi1Pbar = primed
    k2_over_m4 = abs(model.kappaXi1)**2 / model.MXi1**4

    # Eliminating Xi1 in (A.7) gives half the two portal prefactors
    # printed in (D.42), with f_abc = i epsilon_abc / sqrt(2).
    expected_shift = (-portal + primed / np.sqrt(2)) * k2_over_m4
    assert model.alphaOH() - baseline == pytest.approx(expected_shift)
    assert model.alphaOHBox() == pytest.approx(2 * k2_over_m4)
    assert model.alphaOHD() == pytest.approx(4 * k2_over_m4)


@pytest.mark.parametrize("portal", [-0.4, 0.7])
def test_singlet_pure_portal_loop(portal):
    # Expanding the real-scalar determinant with m_S^2(H) = M^2 + 2 lambda_S H†H.
    model = GranadaSMatchingResult()
    model.MS = 1.6
    model.kappaS = model.kappaS3 = model.lam = 0
    model.lambdaS = portal
    loop = 1 / (16 * np.pi**2 * model.MS**2)
    assert model.alphaOH() == pytest.approx(-2 * portal**3 * loop / 3)
    assert model.alphaOHBox() == pytest.approx(-portal**2 * loop / 6)


def test_octet_double_prime_portal_in_yukawa_operators():
    # The bosonic determinant sees |lambda''|². Isolate the terms linear in
    # lambda'' and its conjugate with third-generation exotic Yukawas.
    model = GranadaPhiMatchingResult()
    model.MPhi = 1.7
    third = np.diag([0.0, 0.0, 1.0])
    model.yquPhi = model.yquPhibar = model.ydqPhi = model.ydqPhibar = third
    model.yu = model.yubar = 0.8 * third
    model.yd = model.ydbar = 0.3 * third
    portal = 0.6 + 0.2j

    def coefficients(value):
        model.lambdaHatPrimePrimePhi = value
        model.lambdaHatPrimePrimePhibar = value.conjugate()
        return np.array([model.alphaOuH(2, 2), model.alphaOdH(2, 2)])

    linear = (coefficients(portal) - coefficients(-portal)) / 2
    # The doublet coefficient 1/8 acquires the colour factor 2 C_F = 8/3
    # from Tr(lambda^A lambda^B) and the two quark generators.
    expected = -np.array([0.8 * portal, 0.3 * portal.conjugate()]) / (3 * np.pi**2 * model.MPhi**2)
    np.testing.assert_allclose(linear, expected, rtol=1e-12, atol=1e-15)
