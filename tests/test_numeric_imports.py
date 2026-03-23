import importlib
import inspect
import sys
import types
from pathlib import Path

import pytest


if "wilson" not in sys.modules:
    wilson_stub = types.ModuleType("wilson")
    wilson_stub.wcxf = types.SimpleNamespace()
    sys.modules["wilson"] = wilson_stub


NUMERIC_DIR = Path(__file__).resolve().parents[1] / "lsme" / "numeric"
MODULE_NAMES = sorted(
    f"lsme.numeric.{path.stem}"
    for path in NUMERIC_DIR.glob("*_matching.py")
)


def _matching_classes(module_name):
    module = importlib.import_module(module_name)
    return [
        cls
        for _, cls in inspect.getmembers(module, inspect.isclass)
        if cls.__module__ == module.__name__
    ]


@pytest.mark.parametrize("module_name", MODULE_NAMES)
def test_numeric_matching_classes_import_and_instantiate(module_name):
    classes = _matching_classes(module_name)

    assert classes, f"No classes found in {module_name}"

    for cls in classes:
        instance = cls()
        assert isinstance(instance, cls)


def test_generic_matching_result_import_and_instantiate():
    module = importlib.import_module("lsme.numeric.matchingresult")
    instance = module.GenericMatchingResult(name="TestModel")

    assert instance.name == "TestModel"
