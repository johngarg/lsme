#!/usr/bin/env python3
"""Validate staged classes against black-box evaluations of the original classes.

The reference finite part is reconstructed by polynomial interpolation, without
using MatchMakerParser or the regenerated expressions. Can also write the small
reference fixture used by the ordinary package test suite.
"""

import argparse
import importlib.util
import itertools
import json
from pathlib import Path
import subprocess
import sys
import types

import numpy as np

from regenerate_finite_matching import BASELINE, MODELS, REGULATORS, attributes, methods


def encode(value):
    value = complex(value)
    return [value.real, value.imag]


def load_class(source, label):
    module = types.ModuleType(label)
    exec(compile(source, label, "exec"), module.__dict__)
    name, _ = methods(source)
    return getattr(module, name)


def parameters(model, point):
    rng = np.random.default_rng(20261103 + point)
    # The scale comes first: setting it reloads every SM parameter, so the
    # explicit inputs below must follow it. lam and muH keep the values they
    # had when derived from vev = 0.246 and mH = 0.1251.
    scalars = dict(scale=0.8 if point == 0 else 3.1, g1=0.35 + 0.02 * point,
                   g2=0.63 - 0.01 * point, g3=1.13 + 0.03 * point,
                   lam=0.1251**2 / (2 * 0.246**2), muH=np.sqrt(0.1251**2 / 2))
    arrays = {}
    names = sorted(set(model.exotic_params) | {"yu", "yd", "yl"})
    for offset, name in enumerate(names):
        if name.endswith("bar"):
            continue
        original = getattr(model, name)
        if isinstance(original, np.ndarray):
            value = rng.normal(0, 0.17, original.shape) + 1j * rng.normal(0, 0.11, original.shape)
            arrays[name] = {"real": value.real.tolist(), "imag": value.imag.tolist()}
            if hasattr(model, name + "bar"):
                arrays[name + "bar"] = {"real": value.real.tolist(), "imag": (-value.imag).tolist()}
        else:
            value = 1.7 + 0.4 * point if name.startswith("M") else 0.23 + 0.013 * offset + 0.04 * point
            scalars[name] = value
            if hasattr(model, name + "bar"):
                scalars[name + "bar"] = value
    return {"scalars": scalars, "arrays": arrays}


def configure(model, data):
    for name, value in data["scalars"].items():
        setattr(model, name, value)
    for name, value in data["arrays"].items():
        setattr(model, name, np.asarray(value["real"]) + 1j * np.asarray(value["imag"]))


def evaluate(model, method, indices, epsilon=0.0, pole=0.0, loop=1):
    model.epsilonbar, model.invepsilonbar, model.onelooporder = epsilon, pole, loop
    value = complex(getattr(model, method)(*indices))
    if not np.isfinite(value):
        raise ValueError(f"Nonfinite result: {method}{indices}")
    return value


def reference(model, method, indices):
    samples = {(e, p): evaluate(model, method, indices, e, p)
               for e in (0, 1, 2) for p in (0, 1)}
    base = samples[0, 0]
    pole = samples[0, 1] - base
    mixed1 = samples[1, 1] - samples[1, 0] - pole
    mixed2 = samples[2, 1] - samples[2, 0] - pole
    b1 = 2 * mixed1 - mixed2 / 2
    b2 = mixed2 / 2 - mixed1
    return {
        "finite": base + b1,
        "pole_residue": pole,
        "tree": evaluate(model, method, indices, loop=0),
        "regulator_probe": evaluate(model, method, indices, 0.25, 0) + b1 + 0.25 * b2 + 0.4 * pole,
    }


def actual(model, method, indices):
    finite = evaluate(model, method, indices)
    return {
        "finite": finite,
        "pole_residue": evaluate(model, method, indices, 0, 1) - finite,
        "tree": evaluate(model, method, indices, loop=0),
        "regulator_probe": evaluate(model, method, indices, 0.25, 0.4),
    }


def refresh_varphi_fixture(repo, baseline_class, output):
    """Retain the historical parameters, imposing the new regulator algebra.

    Only the 18 affected entries are updated, using the original class, never
    the candidate. All quartic couplings being tested remain zero here.
    """
    spec = importlib.util.spec_from_file_location("varphi_tests", repo / "tests/test_varphi_matching.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    data = json.loads((repo / "tests/varphi_zero_coupling_baseline.json").read_text())
    source = subprocess.check_output(
        ["git", "-C", str(repo), "show", f"{BASELINE}:lsme/numeric/Granadavarphi_matching.py"], text=True)
    _, functions = methods(source)
    for point in (0, 1):
        model = baseline_class()
        module.configure_baseline_point(model, point)
        model.lambdaHatPrimePrimevarphi = model.lambdaHatPrimePrimevarphibar = 0
        epsilon, pole = model.epsilonbar, model.invepsilonbar
        for name, function in functions.items():
            if not REGULATORS <= attributes(function):
                continue
            indices = (0,) * (len(function.args.args) - 1)
            samples = {(e, p): evaluate(model, name, indices, e, p)
                       for e in (0, 1, 2) for p in (0, 1)}
            b0 = samples[0, 1] - samples[0, 0]
            mixed1 = samples[1, 1] - samples[1, 0] - b0
            mixed2 = samples[2, 1] - samples[2, 0] - b0
            b1, b2 = 2 * mixed1 - mixed2 / 2, mixed2 / 2 - mixed1
            value = evaluate(model, name, indices, epsilon, 0) + b1 + epsilon * b2 + pole * b0
            if abs(value.imag) > 1e-14:
                raise ValueError("Unexpected complex historical fixture")
            data[str(point)][name] = value.real
    output.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--all-flavours", action="store_true")
    parser.add_argument("--write-fixtures", type=Path)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo))
    fixture = {"baseline_commit": BASELINE,
               "prescription": "A(0)+B'(0), from six original-class evaluations; deg_e<=2, deg_p<=1",
               "cases": []}
    report = {"baseline_commit": BASELINE, "all_flavours": args.all_flavours,
              "checks": 0, "max_absolute_error": 0.0, "failures": [], "models": []}
    for model in MODELS:
        filename = model + "_matching.py"
        source = subprocess.check_output(
            ["git", "-C", str(repo), "show", f"{BASELINE}:lsme/numeric/{filename}"], text=True)
        baseline_class = load_class(source, "baseline_" + model)
        candidate_class = load_class((args.candidate_dir / filename).read_text(), "candidate_" + model)
        _, functions = methods(source)
        changed = sorted(name for name, function in functions.items()
                         if REGULATORS <= attributes(function))
        model_checks = 0
        for point in (0, 1):
            old, new = baseline_class(), candidate_class()
            params = parameters(old, point)
            configure(old, params)
            configure(new, params)
            case = {"module": model + "_matching", "class": candidate_class.__name__,
                    "point": point, "parameters": params, "samples": []}
            for name in changed:
                arity = len(functions[name].args.args) - 1
                indices_list = (itertools.product(range(3), repeat=arity) if args.all_flavours
                                else [(0,) * arity, (0, 1, 2, 2)[:arity], (1, 2, 0, 1)[:arity]])
                for indices in indices_list:
                    expected = reference(old, name, indices)
                    observed = actual(new, name, indices)
                    for kind, value in expected.items():
                        error = abs(value - observed[kind])
                        report["max_absolute_error"] = max(report["max_absolute_error"], error)
                        report["checks"] += 1
                        model_checks += 1
                        if not np.isclose(value, observed[kind], rtol=2e-10, atol=2e-12):
                            report["failures"].append(dict(model=model, point=point, method=name,
                                indices=indices, kind=kind, expected=encode(value), actual=encode(observed[kind])))
                    case["samples"].append(dict(method=name, indices=list(indices),
                                                **{kind: encode(value) for kind, value in expected.items()}))
            fixture["cases"].append(case)
        report["models"].append(dict(model=model, methods=len(changed), checks=model_checks))
        print(f"{model}: {model_checks} checks; {len(report['failures'])} cumulative failures", flush=True)
        if args.write_fixtures and model == "Granadavarphi":
            args.write_fixtures.mkdir(parents=True, exist_ok=True)
            refresh_varphi_fixture(repo, baseline_class, args.write_fixtures / "varphi_zero_coupling_baseline.json")
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    if args.write_fixtures:
        (args.write_fixtures / "finite_matching_reference.json").write_text(json.dumps(fixture, indent=2) + "\n")
    if report["failures"]:
        raise SystemExit(f"{len(report['failures'])} validation failures")
    print(f"Passed {report['checks']} comparisons; max absolute error {report['max_absolute_error']:.3g}")


if __name__ == "__main__":
    sys.setrecursionlimit(100000)
    main()
