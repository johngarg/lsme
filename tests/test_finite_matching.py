"""Numerical references from the original classes, independently interpolated.

The fixture uses complex flavor matrices at two mass/scale points. It retains
tree terms, pole residues and an off-limit regulator probe as separate checks.
See tools/validate_finite_matching.py for its reproducible construction.
"""

import importlib
import json
from pathlib import Path

import numpy as np
import pytest


REFERENCE = json.loads(Path(__file__).with_name("finite_matching_reference.json").read_text())


@pytest.mark.parametrize("case", REFERENCE["cases"],
                         ids=lambda case: f"{case['module']}-{case['point']}")
def test_finite_matching_against_original_classes(case):
    module = importlib.import_module("lsme.numeric." + case["module"])
    result = getattr(module, case["class"])()
    for name, value in case["parameters"]["scalars"].items():
        setattr(result, name, value)
    for name, value in case["parameters"]["arrays"].items():
        setattr(result, name, np.asarray(value["real"]) + 1j * np.asarray(value["imag"]))

    for sample in case["samples"]:
        method = getattr(result, sample["method"])
        indices = sample["indices"]
        result.epsilonbar = result.invepsilonbar = 0
        result.onelooporder = 1
        finite = method(*indices)
        result.invepsilonbar = 1
        pole_residue = method(*indices) - finite
        result.invepsilonbar = 0
        result.onelooporder = 0
        tree = method(*indices)
        result.onelooporder = 1
        result.epsilonbar, result.invepsilonbar = 0.25, 0.4
        regulator_probe = method(*indices)
        for name, actual in dict(finite=finite, pole_residue=pole_residue,
                                 tree=tree, regulator_probe=regulator_probe).items():
            expected = complex(*sample[name])
            assert actual == pytest.approx(expected, rel=2e-10, abs=2e-12), (
                case["module"], case["point"], sample["method"], indices, name
            )
