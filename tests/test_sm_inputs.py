"""Physical normalization checks for the fixed MZ Standard Model inputs.

Reference: Huang and Zhou, Phys. Rev. D 103, 016010 (2021), Tables II--V.
https://inspirehep.net/files/71930f37457ebd10742e55c2b36531ea
"""

import numpy as np
import pytest

from lsme.numeric.matchingresult import GEV, GenericMatchingResult


def test_gauge_inputs_reproduce_mz_reference_observables():
    result = GenericMatchingResult("SMInputs")

    # These derived values use the full-SM MZ row of Table V. The weak
    # mass is a tree-level reconstruction, not the measured pole mass.
    weak_mixing_angle = result.g1**2 / (result.g1**2 + result.g2**2)
    weak_mass_gev = result.g2 * result.vev / (2 * GEV)
    alpha_s = result.g3**2 / (4 * np.pi)

    assert weak_mixing_angle == pytest.approx(0.231452935, abs=1e-9)
    assert weak_mass_gev == pytest.approx(80.073, abs=1e-6)
    assert alpha_s == pytest.approx(0.11658642, abs=1e-8)


@pytest.mark.parametrize(
    "coupling,masses_gev",
    [
        ("yu", [0.00123, 0.620, 168.26]),
        ("yd", [0.00267, 0.05316, 2.839]),
        ("yl", [0.00048307, 0.101766, 1.72856]),
    ],
)
def test_yukawas_reproduce_mz_running_masses(coupling, masses_gev):
    result = GenericMatchingResult("SMInputs")
    yukawa = getattr(result, coupling)

    # With <H^0> = v/sqrt(2), the mass matrix must reproduce the full-SM
    # MZ row of Tables II/III, including its diagonal flavor convention.
    reconstructed_masses = yukawa * result.vev / np.sqrt(2)
    np.testing.assert_allclose(
        reconstructed_masses, np.diag(masses_gev) * GEV, rtol=1e-14, atol=0
    )
    np.testing.assert_allclose(
        getattr(result, coupling + "bar"), yukawa.conjugate(), rtol=0, atol=0
    )
