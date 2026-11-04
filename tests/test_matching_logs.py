"""Check matching logarithms against Wilson's one-loop SMEFT beta functions.

For a coefficient i whose tree-level matching vanishes identically,

    dC_i^(1)/d ln mu |_(fixed parameters) = beta_i^EFT[C^(0)],

where C^(0) is the full tree-level coefficient vector, which supplies the
mixing into i. Coefficients with nonzero tree-level matching are not tested:
their identity also needs the UV running of the heavy-field parameters.

The tested coefficients are a selection: tree-level zero both structurally
(every term of the method carries `onelooporder`) and numerically (exactly 0
at several random points), excluding entries far below the largest beta
function.

The redundant ee beta tensor is Fierz-projected before WCxf conversion;
see wilson_beta for the Wilson 2.5.2 caveat.

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
    running = {k: v / (16 * np.pi**2) for k, v in beta(C).items()}
    # Wilson 2.5.2's raw ee beta is symmetric under current exchange but not
    # under Fierz. arrays2wcxf_nonred assumes both. Project the redundant
    # tensor before converting, so the test compares the same Lagrangian.
    ee = running["ee"]
    running["ee"] = (ee + ee.transpose(2, 1, 0, 3)) / 2
    return smeftutil.arrays2wcxf_nonred(running)


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


def mismatches(result, Lambda=None, exclude=(), rtol=1e-9, Gu=None):
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
def test_identical_current_logs_match_wilson_beta(name, operators):
    # Compare the exported logarithms directly; no correction in the test.
    # Count permutation orbits to ensure all multiplicities are exercised.
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
        orbit = {(p, r, s_, t), (s_, t, p, r)}
        if op == "ee":
            orbit |= {(s_, r, p, t), (p, t, s_, r)}
        used.add(len(orbit))
        checked += 1
        assert complex(lhs[k]) == pytest.approx(complex(rhs[k]), rel=1e-8)
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
def test_dipole_logs_match_wilson_beta(name):
    result = getattr(numeric, f"Granada{name}MatchingResult")()
    randomize(result, np.random.default_rng(5))
    checked, bad = mismatches(result)
    assert checked >= 20 and not bad
    lhs = log_derivative(result)
    rhs = wilson_beta(result, tree_level(result))
    targets = [k for k in tree_zero_targets(type(result))
               if k.split("_")[0] in DIPOLE_OPERATORS and abs(rhs.get(k, 0)) > 1e-20]
    assert len(targets) >= 10
    for k in targets:
        assert lhs[k] == pytest.approx(rhs[k], rel=1e-8, abs=1e-22)


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


def test_wilson_ee_beta_projection_preserves_the_operator_sum():
    # With only C_phie^22=1 and g'=1, Wilson's raw beta has
    # beta_1122 = beta_2211 = -1/3, beta_1221 = beta_2112 = 0.
    # All four operators are equal; their total is -2/3, not -4/3.
    result = numeric.GranadaEMatchingResult()
    result.g1, result.g2, result.g3 = 1, 0, 0
    result.yu = result.yd = result.yl = np.zeros((3, 3))
    projected = wilson_beta(result, {"phie_22": 1})
    assert projected["ee_1122"] == pytest.approx(-2 / (3 * 16 * np.pi**2))


def test_omega1_dipole_logs_match_without_diquark_couplings():
    result = numeric.Granadaomega1MatchingResult()
    randomize(result, np.random.default_rng(5))
    for name in ("yqqomega1", "yduomega1"):
        setattr(result, name, np.zeros((3, 3)))
        setattr(result, name + "bar", np.zeros((3, 3)))
    lhs = log_derivative(result)
    rhs = wilson_beta(result, tree_level(result))
    checked = 0
    for key, value in lhs.items():
        if key.split("_")[0] not in DIPOLE_OPERATORS:
            continue
        expected = rhs.get(key, 0)
        if max(abs(value), abs(expected)) > 1e-20:
            checked += 1
        assert value == pytest.approx(expected, rel=1e-8, abs=1e-22)
    assert checked == 36


def test_omega1_diquark_residual_persists_with_symmetric_coupling():
    # Diagnostic of an unresolved matching/running discrepancy, not an
    # expected physical identity. Neither the exporter correction nor the
    # required yqq flavour symmetry fixes the diquark contribution.
    result = numeric.Granadaomega1MatchingResult()
    randomize(result, np.random.default_rng(5))
    result.yqqomega1 = (result.yqqomega1 + result.yqqomega1.T) / 2
    result.yqqomega1bar = result.yqqomega1.conjugate()
    lhs = log_derivative(result)
    rhs = wilson_beta(result, tree_level(result))
    assert abs(rhs["dG_33"]) > 1e-12
    assert lhs["dG_33"] == pytest.approx(-2 * rhs["dG_33"], rel=1e-8, abs=1e-22)
    for key in ("eB_33", "eW_33"):
        assert lhs[key] == pytest.approx(rhs[key], rel=1e-8, abs=1e-22)
