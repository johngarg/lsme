"""Check matching logarithms against Wilson's one-loop SMEFT beta functions.

For a coefficient i whose tree-level matching vanishes identically,

    dC_i^(1)/d ln mu |_(fixed parameters) = beta_i^EFT[C^(0)],

where C^(0) is the full tree-level coefficient vector, which supplies the
mixing into i. Coefficients with nonzero tree-level matching are not tested:
their identity also needs the UV running of the heavy-field parameters.

Both sides are built in Wilson's conventions: C^(0) and the logarithms go
through `to_wcxf`, and the beta functions use the SM parameters the classes
were given. This also tests the Yukawa map, which enters the mixing.
"""

import numpy as np
import pytest
from wilson.run.smeft.beta import beta
from wilson.util import smeftutil

from lsme import numeric
from lsme.numeric.matchingresult import GEV

# to_wcxf exports these with the full-sum normalisation; see
# test_identical_current_operators_need_symmetry_factor.
IDENTICAL_CURRENT_OPERATORS = ("ll", "qq1", "qq3", "uu", "dd", "ee")
LOG_STEP = 0.1


def randomize(result, rng, mass=2.5):
    """Random complex couplings with conjugate barred partners; masses and
    dimensionful kappa couplings real."""
    for name in result.exotic_params:
        if name.endswith("bar") or not hasattr(result, name):
            continue
        shape = np.shape(getattr(result, name))
        if name.startswith("M"):
            setattr(result, name, mass)
            continue
        value = 0.5 * (rng.normal(size=shape) + 1j * rng.normal(size=shape))
        if name.startswith("kappa"):
            value = value.real
        value = value if shape else complex(value)
        setattr(result, name, value)
        if hasattr(result, name + "bar"):
            setattr(result, name + "bar", np.conjugate(value))
    result.scale = mass


def wcxf_coefficients(result):
    return result.to_wcxf().dict


def tree_level(result):
    loop_order = result.onelooporder
    result.onelooporder = 0
    try:
        return wcxf_coefficients(result)
    finally:
        result.onelooporder = loop_order


def log_derivative(result):
    """dC/d ln mu at fixed parameters; one-loop C is linear in ln mu."""
    scale = result._scale
    before = wcxf_coefficients(result)
    try:
        result._scale = scale * np.exp(LOG_STEP)
        after = wcxf_coefficients(result)
    finally:
        result._scale = scale
    return {k: (after[k] - before[k]) / LOG_STEP for k in before}


def wilson_beta(result, tree, Lambda=None):
    """beta_i[C^(0)] in WCxf normalisation, with the classes' SM parameters."""
    C = smeftutil.add_missing(smeftutil.wcxf2arrays_symmetrized(tree))
    C.update({
        "g": result.g2, "gp": result.g1, "gs": result.g3,
        "Lambda": 2 * result.lam if Lambda is None else Lambda,
        "m2": -(result.muH**2).real / GEV**2,
        "Gu": result.yu, "Gd": result.yd, "Ge": result.yl,
    })
    return smeftutil.arrays2wcxf_nonred({k: v / (16 * np.pi**2) for k, v in beta(C).items()})


def identically_zero_at_tree_level(cls, seed, points=3):
    """Coefficients that vanish exactly at tree level at several random points."""
    zero = None
    rng = np.random.default_rng(seed)
    for _ in range(points):
        result = cls()
        randomize(result, rng)
        vanishing = {k for k, v in tree_level(result).items() if v == 0}
        zero = vanishing if zero is None else zero & vanishing
    return zero


def mismatches(result, Lambda=None, exclude=IDENTICAL_CURRENT_OPERATORS, rtol=1e-9):
    targets = identically_zero_at_tree_level(type(result), seed=11)
    lhs = log_derivative(result)
    rhs = wilson_beta(result, tree_level(result), Lambda)
    size = max(abs(complex(rhs.get(k, 0))) for k in targets)
    if size == 0:
        return 0, {}
    bad = {}
    checked = 0
    for k in targets:
        if k.split("_")[0] in exclude:
            continue
        a, b = complex(lhs[k]), complex(rhs.get(k, 0))
        if max(abs(a), abs(b)) < 1e-6 * size:
            continue
        checked += 1
        # Entries far below the largest one lose relative precision to
        # cancellations, so the tolerance has an absolute floor.
        if abs(a - b) > rtol * max(abs(a), abs(b), 1e-3 * size):
            bad[k] = (a, b)
    return checked, bad


# Classes whose logs agree with Wilson's beta for generic complex couplings.
# Not covered: vector models (tree level only), models whose couplings need a
# flavour symmetry (S1, S2, Upsilon, Xi1, Omega2, Omega4, omega2, omega4,
# zeta), and the cases pinned down by the tests below (S, N, Sigma, dipoles).
AGREEING_CLASSES = [
    "D", "Delta1", "Delta3", "E", "ND", "Pi1", "Q1", "Q5", "Q7",
    "Sigma1", "SigmaD", "T1", "T2", "U", "Xi",
]


@pytest.mark.parametrize("name", AGREEING_CLASSES)
def test_matching_logs_match_wilson_beta(name):
    cls = getattr(numeric, f"Granada{name}MatchingResult")
    result = cls()
    randomize(result, np.random.default_rng(5))
    checked, bad = mismatches(result)
    assert checked >= 20
    assert not bad


def test_singlet_scalar_logs_use_uv_quartic():
    # The singlet's trilinear shifts the EFT quartic at tree level,
    # lambda_EFT = lam − kappa^2/(2 M^2). The class's lam is the UV quartic:
    # its logs agree with Wilson's beta only after this shift, so feeding it
    # Wilson's (EFT) quartic is wrong at O(kappa^2/M^2).
    result = numeric.GranadaSMatchingResult()
    randomize(result, np.random.default_rng(7))
    shift = result.kappaS**2 / (2 * result.MS**2)
    checked, bad = mismatches(result, Lambda=2 * (result.lam - shift))
    assert checked >= 20 and not bad
    _, unshifted = mismatches(result)
    assert unshifted


@pytest.mark.xfail(strict=True, reason=(
    "to_wcxf exports identical-current four-fermion coefficients with the "
    "full-sum normalisation; WCxf lists each distinct operator once, so "
    "entries with (pr) != (st) need a factor 2"))
def test_identical_current_operators_need_symmetry_factor():
    result = numeric.GranadaEMatchingResult()
    randomize(result, np.random.default_rng(5))
    _, bad = mismatches(result, exclude=())
    assert not bad


@pytest.mark.parametrize("wrong_map", [
    lambda y: y.T,
    lambda y: np.conjugate(y),
])
def test_wrong_up_yukawa_map_breaks_logs(wrong_map):
    result = numeric.GranadaUMatchingResult()
    randomize(result, np.random.default_rng(5))
    result.yu = wrong_map(result.yu)
    result.yubar = np.conjugate(result.yu)
    lhs = log_derivative(result)
    # Wilson's beta keeps the correct Gu, as the table supplies it.
    from lsme.numeric.sm_running import sm_parameters
    correct = sm_parameters(result.scale)
    reference = numeric.GranadaUMatchingResult()
    reference.scale = result.scale
    reference.yu = correct["yu"]
    rhs = wilson_beta(reference, tree_level(result))
    assert any(
        abs(complex(lhs[k]) - complex(rhs.get(k, 0))) > 1e-3 * abs(complex(rhs.get(k, 0)))
        for k in lhs
        if k.startswith("uphi_") and abs(complex(rhs.get(k, 0))) > 0
    )


@pytest.mark.xfail(strict=True, reason=(
    "dipole coefficients are exported with the opposite sign to Wilson's "
    "convention: their logs equal exactly minus Wilson's beta, consistent with "
    "opposite signs of the gauge couplings in the covariant derivative"))
@pytest.mark.parametrize("name", ["Phi", "Pi7", "varphi"])
def test_dipole_logs_match_wilson_beta(name):
    result = getattr(numeric, f"Granada{name}MatchingResult")()
    randomize(result, np.random.default_rng(5))
    _, bad = mismatches(result)
    assert not bad


def test_dipole_logs_are_minus_wilson_beta():
    result = numeric.GranadaPhiMatchingResult()
    randomize(result, np.random.default_rng(5))
    _, bad = mismatches(result)
    dipoles = ("dB", "dG", "dW", "uB", "uG", "uW", "eB", "eW")
    assert bad and all(k.split("_")[0] in dipoles for k in bad)
    for a, b in bad.values():
        assert a == pytest.approx(-b, rel=1e-9)


@pytest.mark.parametrize("name", ["N", "Sigma"])
def test_majorana_logs_exceed_wilson_beta_at_quartic_order(name):
    # For the Majorana seesaw fields the logs carry extra terms quartic in
    # the heavy coupling, i.e. quadratic in the dimension-five coefficient,
    # which Wilson's dimension-six beta functions do not include.
    def discrepancy(factor):
        result = getattr(numeric, f"Granada{name}MatchingResult")()
        randomize(result, np.random.default_rng(5))
        coupling = f"lambda{name}"
        setattr(result, coupling, factor * getattr(result, coupling))
        setattr(result, coupling + "bar", np.conjugate(getattr(result, coupling)))
        _, bad = mismatches(result)
        return max(abs(a - b) for a, b in bad.values())

    assert discrepancy(0.5) / discrepancy(0.25) == pytest.approx(16, rel=0.05)
