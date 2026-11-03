"""Check the precomputed table of Wilson's pure-SM parameters."""

import json
import time
import warnings

import numpy as np
import pytest
from numpy.polynomial import chebyshev

from lsme.numeric import sm_table
from lsme.numeric.sm_running import loaded_table, wilson_sm_lagrangian
from lsme.numeric.sm_table import (
    MZ_GEV,
    MU_MAX_GEV,
    SMTable,
    VALIDATION_RTOL,
    live_wilson_sm_lagrangian,
    load_sm_table,
    physical_quantities,
    wilson_fingerprint,
)

OFF_GRID_SCALES = [MZ_GEV, 137.0, 512.3, 2.9e3, 1.7e4, 8.8e4, 3.1e5, MU_MAX_GEV]


def _relative_error(params, reference):
    approx = physical_quantities(params)
    exact = physical_quantities(reference)
    return np.max(np.abs(approx - exact) / np.abs(exact))


def _stale_copy(tmp_path, **metadata):
    table = loaded_table()
    stale = SMTable(table.coefficients, table.mu_max_gev, {**table.metadata, **metadata})
    path = tmp_path / "sm_table.npz"
    stale.save(path)
    return path


@pytest.fixture(scope="module")
def rebuilt_table(tmp_path_factory):
    path = _stale_copy(tmp_path_factory.mktemp("stale"), wilson_version="0.0.0")
    with pytest.warns(RuntimeWarning, match="different wilson_version"):
        return load_sm_table(path)


def test_shipped_table_loads_and_matches_installed_wilson():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        table = load_sm_table()
    for key, value in wilson_fingerprint().items():
        assert table.metadata[key] == value
    assert table.metadata["mu_range_gev"] == [MZ_GEV, MU_MAX_GEV]
    assert table.mu_max_gev == MU_MAX_GEV


@pytest.mark.parametrize("scale", OFF_GRID_SCALES)
def test_table_agrees_with_live_wilson(scale):
    # Couplings, each Yukawa singular value relative to itself, CKM
    # magnitudes and the Jarlskog invariant.
    assert _relative_error(wilson_sm_lagrangian(scale), live_wilson_sm_lagrangian(scale)) < VALIDATION_RTOL


def test_gauge_running_matches_one_loop_beta_functions():
    table = loaded_table()
    derivative = chebyshev.chebder(table.coefficients)
    log_range = np.log(table.mu_max_gev / MZ_GEV)
    b = {"gp": 41 / 6, "g": -19 / 6, "gs": -7}
    for scale in np.geomspace(2 * MZ_GEV, table.mu_max_gev / 2, 9):
        params = wilson_sm_lagrangian(scale)
        dvalues = chebyshev.chebval(table._x(scale), derivative) * 2 / log_range
        for n, key in enumerate(("g", "gp", "gs")):
            expected = b[key] * params[key]**3 / (16 * np.pi**2)
            assert dvalues[n] == pytest.approx(expected, rel=1e-2)


def test_flavour_entries_vary_smoothly():
    scales = np.geomspace(MZ_GEV, MU_MAX_GEV, 400)
    for key in ("Gu", "Gd", "Ge"):
        values = np.array([wilson_sm_lagrangian(s)[key] for s in scales])
        size = np.abs(values).max()
        jumps = np.abs(np.diff(values, axis=0)).max()
        assert jumps < 1e-2 * size


@pytest.mark.parametrize("scale", [np.nan, np.inf, MZ_GEV * (1 - 1e-12), 0.0, -1.0, MU_MAX_GEV * (1 + 1e-12)])
def test_scales_outside_table_raise(scale):
    with pytest.raises(ValueError):
        wilson_sm_lagrangian(scale)


def test_endpoints_are_accepted():
    for scale in (MZ_GEV, MU_MAX_GEV):
        assert np.isfinite(wilson_sm_lagrangian(scale)["gs"])


def test_stale_wilson_version_rebuilds_in_memory(rebuilt_table):
    for key, value in wilson_fingerprint().items():
        assert rebuilt_table.metadata[key] == value
    scale = 4.2e3
    assert _relative_error(rebuilt_table(scale), live_wilson_sm_lagrangian(scale)) < VALIDATION_RTOL


@pytest.mark.parametrize("reason,make_path", [
    ("different smpar_sha256", lambda tmp: _stale_copy(tmp, smpar_sha256="0" * 64)),
    ("could not be read", lambda tmp: tmp / "missing.npz"),
])
def test_unusable_table_warns_and_rebuilds(monkeypatch, tmp_path, reason, make_path):
    sentinel = loaded_table()
    monkeypatch.setattr(sm_table, "build_sm_table", lambda: sentinel)
    path = make_path(tmp_path)
    with pytest.warns(RuntimeWarning, match=reason):
        assert load_sm_table(path) is sentinel
    assert not (tmp_path / "missing.npz").exists()


def test_regenerated_table_reproduces_shipped_table(rebuilt_table):
    for scale in np.geomspace(MZ_GEV, MU_MAX_GEV, 25):
        assert _relative_error(rebuilt_table(scale), wilson_sm_lagrangian(scale)) < VALIDATION_RTOL


def test_lookup_and_load_are_fast():
    times = []
    for scale in np.geomspace(MZ_GEV, MU_MAX_GEV, 1000):
        start = time.perf_counter()
        wilson_sm_lagrangian(scale)
        times.append(time.perf_counter() - start)
    assert np.median(times) < 100e-6

    start = time.perf_counter()
    load_sm_table()
    assert time.perf_counter() - start < 50e-3


def test_lookups_return_independent_values():
    first = wilson_sm_lagrangian(1e3)
    first["Gu"][2, 2] = 0
    first["gs"] = 0
    second = wilson_sm_lagrangian(1e3)
    assert second["Gu"][2, 2] != 0
    assert second["gs"] != 0
    assert all(isinstance(second[key], float) for key in ("g", "gp", "gs", "Lambda", "m2"))
    assert all(second[key].shape == (3, 3) for key in ("Gu", "Gd", "Ge"))


def test_saved_metadata_is_plain_json(tmp_path):
    path = tmp_path / "table.npz"
    loaded_table().save(path)
    with np.load(path, allow_pickle=False) as data:
        assert json.loads(str(data["metadata"]))["degree"] == loaded_table().metadata["degree"]
