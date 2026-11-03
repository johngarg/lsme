#!/usr/bin/env python3
"""Wilson's pure-SM MS-bar parameters at a given scale, from the shipped table.

The table is loaded when this module is imported; see `sm_table` for its
contents, conventions and accuracy.

`sm_parameters` maps Wilson's parameters to the LSME (MatchMaker) SM
conventions and TeV units:

    LSME:   L ⊃ −muH^2 H†H − lam (H†H)^2
                − (yd Q̄ D H + yl L̄ E H + yu Q̄ U H̃ + h.c.)
    Wilson: L ⊃ m2 φ†φ − (Lambda/2)(φ†φ)^2
                − (Gd q̄ d φ + Ge l̄ e φ + Gu q̄ u φ̃ + h.c.)

    | LSME                | Wilson                          |
    |---------------------|---------------------------------|
    | g1, g2, g3          | gp, g, gs                       |
    | lam                 | Lambda / 2                      |
    | muH**2              | −m2 (×1e-6, GeV^2 → TeV^2)      |
    | yu, yd, yl          | Gu, Gd, Ge (same index order)   |
    | yubar, ydbar, ylbar | conj(Gu), conj(Gd), conj(Ge)    |

Wilson's m2 is positive, so muH**2 is negative and `muH` is imaginary,
muH = 1j * sqrt(m2). The matching classes only use (self.muH)**(2), which is
then exactly real and negative. That does not make coefficients real in
general: they are complex whenever other inputs are, e.g. through the CKM
phase in yu. The Higgs-mass and Yukawa maps agree with the
MatchMakerEFT-to-Matchete conversion in Matchete's E_VLL example
(mu^2 → −mu^2, lambda → lambda/2, Yukawas unchanged), whose Matchete SM
Lagrangian has Wilson's form. Gu carries the CKM matrix in Wilson's
down-aligned basis.
"""

import numpy as np

from .sm_table import load_sm_table

TEV_IN_GEV = 1e3
SM_PARAMETER_NAMES = ("g1", "g2", "g3", "lam", "muH", "yu", "yd", "yl", "yubar", "ydbar", "ylbar")

_TABLE = load_sm_table()


def wilson_sm_lagrangian(scale_gev):
    """Return Wilson's pure-SM parameters at `scale_gev` from the table.

    Keys g, gp, gs, Lambda, m2 (floats; m2 in GeV**2) and Gu, Gd, Ge (complex
    3x3 arrays), in Wilson's conventions. Valid for MZ <= scale_gev <= the
    table's upper limit (1000 TeV by default); other scales raise ValueError.
    """
    return _TABLE(scale_gev)


def sm_parameters(scale_tev):
    """Return the LSME SM parameters at `scale_tev` (in TeV), as fresh objects.

    Keys are SM_PARAMETER_NAMES; see the module docstring for the map from
    Wilson's conventions. Scales outside the table raise ValueError.
    """
    wilson = wilson_sm_lagrangian(float(scale_tev) * TEV_IN_GEV)
    muH2 = -wilson["m2"] / TEV_IN_GEV**2
    yukawas = {"yu": wilson["Gu"], "yd": wilson["Gd"], "yl": wilson["Ge"]}
    return {
        "g1": wilson["gp"],
        "g2": wilson["g"],
        "g3": wilson["gs"],
        "lam": wilson["Lambda"] / 2,
        "muH": float(np.sqrt(muH2)) if muH2 >= 0 else complex(0.0, np.sqrt(-muH2)),
        **yukawas,
        **{name + "bar": np.conjugate(y) for name, y in yukawas.items()},
    }


def loaded_table():
    """Return the SMTable loaded at import."""
    return _TABLE
