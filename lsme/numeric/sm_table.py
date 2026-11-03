#!/usr/bin/env python3
"""Build, save and load the table of Wilson's pure-SM MS-bar parameters.

The parameters are those Wilson extracts for an empty set of SMEFT
coefficients, `SMEFT(wcxf.WC("SMEFT", "Warsaw", mu, {})).C_in`, in Wilson's
own conventions and GeV units:

    L ⊃ m2 φ†φ − (Lambda/2)(φ†φ)^2 − (q̄ Gu u φ̃ + q̄ Gd d φ + l̄ Ge e φ + h.c.),

with gauge couplings g (SU(2)), gp (U(1)_Y) and gs (SU(3)). Gu, Gd and Ge are
in Wilson's default basis: Gd and Ge diagonal, Gu carrying the CKM matrix in
the standard phase convention.

A live extraction costs ~25 ms per scale and carries Wilson's per-scale
integration noise (rtol 1e-3). Instead, a least-squares Chebyshev fit in
t = ln(mu/MZ) is precomputed by `scripts/build_sm_table.py`, shipped as
`data/sm_table.npz`, and loaded by `lsme.numeric.sm_running` at import. A
lookup costs tens of microseconds and varies smoothly with the scale. Its accuracy is
limited by Wilson's integration tolerance, not by the fit.

No mapping to LSME conventions or TeV units is applied here.
"""

import hashlib
import json
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from numpy.polynomial import chebyshev

MZ_GEV = 91.1876
MU_MAX_GEV = 1e6
SCALAR_KEYS = ("g", "gp", "gs", "Lambda", "m2")
MATRIX_KEYS = ("Gu", "Gd", "Ge")
COMPONENTS = SCALAR_KEYS + tuple(
    f"{key}_{part}_{i}{j}"
    for key in MATRIX_KEYS
    for part in ("re", "im")
    for i in range(3)
    for j in range(3)
)
TABLE_PATH = Path(__file__).with_name("data") / "sm_table.npz"
VALIDATION_RTOL = 1e-3
SCALAR_IMAG_RTOL = 1e-10


def _check_scale(scale_gev, mu_max_gev=None):
    scale_gev = float(scale_gev)
    if not np.isfinite(scale_gev):
        raise ValueError(f"Scale must be finite, got {scale_gev}")
    if scale_gev < MZ_GEV or (mu_max_gev is not None and scale_gev > mu_max_gev):
        upper = "" if mu_max_gev is None else f", {mu_max_gev:g}"
        raise ValueError(f"Scale {scale_gev:g} GeV is outside [{MZ_GEV}{upper}] GeV")
    return scale_gev


def wilson_fingerprint():
    """Return the Wilson version and a hash of its default SM inputs."""
    import wilson
    from wilson.run.smeft import smpar

    inputs = repr(sorted(smpar.p.items())).encode()
    return {
        "wilson_version": wilson.__version__,
        "smpar_sha256": hashlib.sha256(inputs).hexdigest(),
    }


def live_wilson_sm_lagrangian(scale_gev):
    """Extract Wilson's pure-SM parameters at `scale_gev` (slow, ~25 ms).

    Used to build and check the table, not for repeated evaluation.
    """
    from wilson import wcxf
    from wilson.run.smeft import SMEFT

    scale_gev = _check_scale(scale_gev)
    C = SMEFT(wcxf.WC("SMEFT", "Warsaw", scale_gev, {})).C_in
    params = {}
    for key in SCALAR_KEYS:
        value = complex(C[key])
        if abs(value.imag) > SCALAR_IMAG_RTOL * abs(value.real):
            raise ValueError(f"Wilson returned complex {key} = {value} at {scale_gev:g} GeV")
        params[key] = value.real
    for key in MATRIX_KEYS:
        params[key] = np.array(C[key], dtype=complex)
    return params


def _flatten(params):
    matrices = [np.asarray(params[key], dtype=complex) for key in MATRIX_KEYS]
    return np.concatenate(
        [[params[key] for key in SCALAR_KEYS]]
        + [part for M in matrices for part in (M.real.ravel(), M.imag.ravel())]
    )


def _unflatten(values):
    params = {key: float(values[n]) for n, key in enumerate(SCALAR_KEYS)}
    offset = len(SCALAR_KEYS)
    for key in MATRIX_KEYS:
        real = values[offset:offset + 9]
        imag = values[offset + 9:offset + 18]
        params[key] = (real + 1j * imag).reshape(3, 3)
        offset += 18
    return params


def physical_quantities(params):
    """Return the quantities used to judge the table's accuracy.

    These are the five scalars, the nine Yukawa singular values (each judged
    relative to itself, so light-fermion Yukawas are not swamped by the top),
    the CKM magnitudes, and the Jarlskog invariant.
    """
    Uu, yu, _ = np.linalg.svd(params["Gu"])
    Ud, yd, _ = np.linalg.svd(params["Gd"])
    ye = np.linalg.svd(params["Ge"], compute_uv=False)
    # Order generations by increasing mass.
    V = Uu[:, ::-1].conj().T @ Ud[:, ::-1]
    jarlskog = np.imag(V[0, 1] * V[1, 2] * np.conj(V[0, 2]) * np.conj(V[1, 1]))
    return np.concatenate([
        [params[key] for key in SCALAR_KEYS],
        yu[::-1], yd[::-1], ye[::-1],
        np.abs(V).ravel(),
        [jarlskog],
    ])


def _max_relative_error(params, reference):
    approx = physical_quantities(params)
    exact = physical_quantities(reference)
    return float(np.max(np.abs(approx - exact) / np.abs(exact)))


class SMTable:
    """Chebyshev representation of Wilson's pure-SM parameters in ln(mu)."""

    def __init__(self, coefficients, mu_max_gev, metadata):
        self.coefficients = np.asarray(coefficients, dtype=float)
        self.mu_max_gev = float(mu_max_gev)
        self.metadata = dict(metadata)
        self._log_range = np.log(self.mu_max_gev / MZ_GEV)
        if self.coefficients.shape[1:] != (len(COMPONENTS),):
            raise ValueError(f"Unexpected coefficient shape {self.coefficients.shape}")

    def _x(self, scale_gev):
        return 2 * np.log(scale_gev / MZ_GEV) / self._log_range - 1

    def __call__(self, scale_gev):
        """Return Wilson's raw SM parameters at `scale_gev`, as fresh objects."""
        scale_gev = _check_scale(scale_gev, self.mu_max_gev)
        return _unflatten(chebyshev.chebval(self._x(scale_gev), self.coefficients))

    def save(self, path=TABLE_PATH):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            path,
            coefficients=self.coefficients,
            components=np.array(COMPONENTS),
            mu_range_gev=np.array([MZ_GEV, self.mu_max_gev]),
            metadata=np.array(json.dumps(self.metadata, sort_keys=True)),
        )

    @classmethod
    def load(cls, path=TABLE_PATH):
        with np.load(path, allow_pickle=False) as data:
            if tuple(data["components"]) != COMPONENTS:
                raise ValueError(f"Component layout in {path} does not match this module")
            mu_min_gev, mu_max_gev = data["mu_range_gev"]
            if mu_min_gev != MZ_GEV:
                raise ValueError(f"Table in {path} starts at {mu_min_gev} GeV, not MZ")
            return cls(data["coefficients"], mu_max_gev, json.loads(str(data["metadata"])))


def _chebyshev_nodes(n):
    return np.cos(np.pi * (np.arange(n) + 0.5) / n)[::-1]


def _check_smooth(values):
    """Reject a grid whose entries jump between neighbouring nodes, e.g. from
    a change of flavour phase convention. Entries that vanish by the basis
    choice are judged against their matrix, not against their rounding noise.
    """
    scale = np.abs(values).max(axis=0)
    offset = len(SCALAR_KEYS)
    for _ in MATRIX_KEYS:
        block = slice(offset, offset + 18)
        scale[block] = np.maximum(scale[block], 1e-12 * np.abs(values[:, block]).max())
        offset += 18
    jumps = np.abs(np.diff(values, 2, axis=0)) / scale
    worst = int(np.unravel_index(jumps.argmax(), jumps.shape)[1])
    if jumps.max() > 0.05:
        raise ValueError(
            f"Wilson's SM parameters are not smooth on the grid "
            f"(component {COMPONENTS[worst]}, relative jump {jumps.max():.2g})"
        )
    return float(jumps.max())


def build_sm_table(mu_max_gev=MU_MAX_GEV, n_points=64, n_validation=20,
                   degrees=range(4, 21), seed=0):
    """Fit Wilson's pure-SM parameters on [MZ, mu_max_gev] and validate them.

    The degree is the smallest whose held-out error is within a factor two of
    the best candidate; the held-out reference itself carries Wilson's
    integration noise. Raises ValueError if the grid is not smooth or any
    physical quantity misses VALIDATION_RTOL at a held-out scale.
    """
    mu_max_gev = _check_scale(mu_max_gev)
    log_range = np.log(mu_max_gev / MZ_GEV)
    x_nodes = _chebyshev_nodes(n_points)
    nodes = [live_wilson_sm_lagrangian(MZ_GEV * np.exp((x + 1) * log_range / 2)) for x in x_nodes]
    values = np.array([_flatten(params) for params in nodes])
    max_jump = _check_smooth(values)

    x_held = np.random.default_rng(seed).uniform(-1, 1, n_validation)
    held_scales = MZ_GEV * np.exp((x_held + 1) * log_range / 2)
    held = [live_wilson_sm_lagrangian(scale) for scale in held_scales]

    errors = {}
    fits = {}
    for degree in degrees:
        fits[degree] = chebyshev.chebfit(x_nodes, values, degree)
        errors[degree] = max(
            _max_relative_error(_unflatten(chebyshev.chebval(x, fits[degree])), reference)
            for x, reference in zip(x_held, held)
        )
    best = min(errors.values())
    degree = min(d for d, error in errors.items() if error <= 2 * best)
    if errors[degree] > VALIDATION_RTOL:
        raise ValueError(
            f"Table misses the validation tolerance {VALIDATION_RTOL:g}: "
            f"held-out error {errors[degree]:.2g} at degree {degree}"
        )

    metadata = {
        **wilson_fingerprint(),
        "numpy_version": np.__version__,
        "mu_range_gev": [MZ_GEV, mu_max_gev],
        "n_points": n_points,
        "degree": degree,
        "n_validation": n_validation,
        "max_heldout_relative_error": errors[degree],
        "max_grid_relative_jump": max_jump,
        "phase_convention": "Wilson default basis (ckmutil standard CKM phases)",
        "built": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    return SMTable(fits[degree], mu_max_gev, metadata)


def load_sm_table(path=TABLE_PATH):
    """Load the shipped table, rebuilding it in memory if it is stale.

    A missing or unreadable file, or a different Wilson version or default
    SM inputs, triggers one warning and an in-memory rebuild (~2 s). Nothing
    is written to disk here.
    """
    try:
        table = SMTable.load(path)
    except (OSError, ValueError, KeyError) as exc:
        reason = f"could not be read ({exc})"
    else:
        current = wilson_fingerprint()
        stale = {key for key, value in current.items() if table.metadata.get(key) != value}
        if not stale:
            return table
        reason = "was built for a different " + " and ".join(sorted(stale))
    warnings.warn(
        f"SM parameter table {path} {reason}; rebuilding it in memory. "
        "Regenerate it with scripts/build_sm_table.py.",
        RuntimeWarning,
        stacklevel=2,
    )
    return build_sm_table()
