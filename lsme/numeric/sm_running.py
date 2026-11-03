#!/usr/bin/env python3
"""Wilson's pure-SM MS-bar parameters at a given scale, from the shipped table.

The table is loaded when this module is imported; see `sm_table` for its
contents, conventions and accuracy.
"""

from .sm_table import load_sm_table

_TABLE = load_sm_table()


def wilson_sm_lagrangian(scale_gev):
    """Return Wilson's pure-SM parameters at `scale_gev` from the table.

    Keys g, gp, gs, Lambda, m2 (floats; m2 in GeV**2) and Gu, Gd, Ge (complex
    3x3 arrays), in Wilson's conventions. Valid for MZ <= scale_gev <= the
    table's upper limit (1000 TeV by default); other scales raise ValueError.
    """
    return _TABLE(scale_gev)


def loaded_table():
    """Return the SMTable loaded at import."""
    return _TABLE
