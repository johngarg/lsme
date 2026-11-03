"""Check matching logarithms against Wilson's one-loop SMEFT beta functions.

For a coefficient i whose tree-level matching vanishes identically,

    dC_i^(1)/d ln mu |_(fixed parameters) = beta_i^EFT[C^(0)],

where C^(0) is the full tree-level coefficient vector, which supplies the
mixing into i. Coefficients with nonzero tree-level matching are not tested:
their identity also needs the UV running of the heavy-field parameters.

The tested coefficients are a selection: tree-level zero both structurally
(every term of the method carries `onelooporder`) and numerically (exactly 0
at several random points), excluding the identical-current four-fermion
families and entries far below the largest beta function.

Both sides are built in Wilson's conventions: C^(0) and the logarithms go
through `to_wcxf`, and the beta functions use the SM parameters the classes
were given. This also tests the Yukawa map, which enters the mixing.
"""

import ast
import functools
import inspect

import numpy as np
import pytest
from wilson.run.smeft.beta import beta
from wilson.util import smeftutil

from lsme import numeric
from lsme.numeric.matchingresult import GEV

# to_wcxf exports these with MatchMaker's symmetric full-sum normalisation,
# not WCxf's; see test_identical_current_logs_match_after_symmetry_factors.
IDENTICAL_CURRENT_OPERATORS = ("ll", "qq1", "qq3", "uu", "dd", "ee")
DIPOLE_OPERATORS = ("dB", "dG", "dW", "uB", "uG", "uW", "eB", "eW")
LOG_STEP = 0.1


def randomize(result, rng, mass=2.5):
    """Random complex couplings with conjugate barred partners; masses and
    dimensionful kappa couplings real."""
    missing = [name for name in result.exotic_params if not hasattr(result, name)]
    if missing:
        raise AttributeError(f"{type(result).__name__} lists missing parameters {missing}")
    for name in result.exotic_params:
        if name.endswith("bar"):
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


def wilson_beta(result, tree, Lambda=None, Gu=None):
    """beta_i[C^(0)] in WCxf normalisation, with the classes' SM parameters
    unless Lambda or Gu are given."""
    C = smeftutil.add_missing(smeftutil.wcxf2arrays_symmetrized(tree))
    C.update({
        "g": result.g2, "gp": result.g1, "gs": result.g3,
        "Lambda": 2 * result.lam if Lambda is None else Lambda,
        "m2": -(result.muH**2).real / GEV**2,
        "Gu": result.yu if Gu is None else Gu, "Gd": result.yd, "Ge": result.yl,
    })
    return smeftutil.arrays2wcxf_nonred({k: v / (16 * np.pi**2) for k, v in beta(C).items()})


def _operands(node, ops):
    """Operands of a chain of the given binary operators (iterative: the
    generated expressions are too long to walk recursively)."""
    out, stack = [], [node]
    while stack:
        n = stack.pop()
        if isinstance(n, ast.BinOp) and isinstance(n.op, ops):
            stack += [n.right, n.left]
        elif isinstance(n, ast.UnaryOp) and isinstance(n.op, (ast.USub, ast.UAdd)):
            stack.append(n.operand)
        else:
            out.append(n)
    return out


def _vanishes_at_tree_level(node):
    if isinstance(node, ast.Constant) and node.value == 0:
        return True
    for term in _operands(node, (ast.Add, ast.Sub)):
        factors = _operands(term, (ast.Mult, ast.Div))
        if any(isinstance(f, ast.Attribute) and f.attr == "onelooporder" for f in factors):
            continue
        if any(_operands(f, (ast.Add, ast.Sub)) != [f] and _vanishes_at_tree_level(f) for f in factors):
            continue
        return False
    return True


@functools.lru_cache(maxsize=None)
def structurally_zero_at_tree_level(cls):
    """Coefficients whose method has every term multiplied by onelooporder."""
    module = ast.parse(inspect.getsource(inspect.getmodule(cls)))
    class_node = next(n for n in module.body if isinstance(n, ast.ClassDef) and n.name == cls.__name__)
    returns = {}
    for node in class_node.body:
        if isinstance(node, ast.FunctionDef):
            found = [n for n in ast.walk(node) if isinstance(n, ast.Return)]
            returns[node.name] = found[0].value if len(found) == 1 else None
    zero = set()
    for op in cls().coeff_dict():
        value = returns.get(cls._warsaw_method(op)[0])
        if value is not None and _vanishes_at_tree_level(value):
            zero.add(op)
    return frozenset(zero)


@functools.lru_cache(maxsize=None)
def numerically_zero_at_tree_level(cls, seed, points=3):
    """Coefficients that vanish exactly at tree level at several random points."""
    zero = None
    rng = np.random.default_rng(seed)
    for _ in range(points):
        result = cls()
        randomize(result, rng)
        vanishing = {k for k, v in tree_level(result).items() if v == 0}
        zero = vanishing if zero is None else zero & vanishing
    return frozenset(zero)


def tree_zero_targets(cls):
    structural = structurally_zero_at_tree_level(cls)
    numerical = numerically_zero_at_tree_level(cls, seed=11)
    assert structural <= numerical, sorted(structural - numerical)[:5]
    return structural


def mismatches(result, Lambda=None, exclude=IDENTICAL_CURRENT_OPERATORS, rtol=1e-9, Gu=None):
    targets = tree_zero_targets(type(result))
    lhs = log_derivative(result)
    rhs = wilson_beta(result, tree_level(result), Lambda, Gu)
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


@pytest.mark.parametrize("name,operators", [
    ("E", ("ll",)), ("D", ("qq1", "qq3")), ("Q1", ("uu", "dd")),
    ("Delta1", ("ee",)), ("Delta3", ("ee",)),
])
def test_identical_current_logs_match_after_symmetry_factors(name, operators):
    # Multiplying the exported values by Wilson's own WCxf symmetry factors
    # (1 or 2 for ll, qq1, qq3, uu, dd; 1, 2 or 4 for ee, which has an extra
    # Fierz symmetry) restores agreement: the conversion the exporter lacks.
    factors = {"ee": smeftutil._d_6}
    result = getattr(numeric, f"Granada{name}MatchingResult")()
    randomize(result, np.random.default_rng(5))
    targets = tree_zero_targets(type(result))
    lhs = log_derivative(result)
    rhs = wilson_beta(result, tree_level(result))
    used = set()
    checked = 0
    for k in targets:
        op, _, indices = k.partition("_")
        if op not in operators or abs(complex(rhs.get(k, 0))) < 1e-20:
            continue
        p, r, s_, t = (int(c) - 1 for c in indices)
        factor = factors.get(op, smeftutil._d_4)[p, r, s_, t]
        used.add(float(factor))
        checked += 1
        assert complex(lhs[k]) * factor == pytest.approx(complex(rhs[k]), rel=1e-8)
    assert checked >= 10
    assert used >= ({1.0, 2.0, 4.0} if "ee" in operators else {1.0, 2.0})


@pytest.mark.parametrize("name,wrong_map,breaks", [
    ("U", lambda y: y, False),
    ("U", lambda y: y.T, True),
    ("U", lambda y: np.conjugate(y), True),
])
def test_up_yukawa_map_is_checked_by_logs(name, wrong_map, breaks):
    # The class uses the mapped Yukawa; Wilson's beta keeps the table's Gu.
    # Only tree-level-zero coefficients are compared, so the identity map is a
    # control that must agree.
    result = getattr(numeric, f"Granada{name}MatchingResult")()
    randomize(result, np.random.default_rng(5))
    Gu = result.yu.copy()
    result.yu = wrong_map(Gu)
    result.yubar = np.conjugate(result.yu)
    checked, bad = mismatches(result, Gu=Gu)
    assert checked >= 20
    assert bool(bad) == breaks


@pytest.mark.parametrize("name", ["Phi", "Pi7", "varphi"])
def test_dipole_logs_are_minus_wilson_beta(name):
    result = getattr(numeric, f"Granada{name}MatchingResult")()
    randomize(result, np.random.default_rng(5))
    _, bad = mismatches(result)
    assert bad and all(k.split("_")[0] in DIPOLE_OPERATORS for k in bad)
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
