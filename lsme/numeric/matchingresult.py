#!/usr/bin/env python3

from dataclasses import dataclass
from typing import List
import numpy as np

from wilson import wcxf

GEV = 1e-3
MZ = 91.1876 * GEV
BARYON_VIOLATING_WARSAW_OPERATORS = frozenset({"duql", "duue", "qqql", "qque"})

class GenericMatchingResult:
    """Class for storing MatchMaker matching results. Model parameters, direct
    coefficient methods and coeff_dict use TeV units; to_wcxf exports in GeV
    units.
    Scale is initialised to 1 TeV.
    Gauge and Yukawa defaults are fixed full-SM MS-bar inputs at MZ;
    changing scale does not evolve these inputs.

    """

    def __init__(self, name: str, scale: float = 1.0):
        self.name = name
        self.scale = scale
        self.exotic_params = {}
        self.sm_params = ['g1', 'g2', 'g3', 'mH', 'mZ', 'vev', 'yd', 'ydbar', 'yl', 'ylbar', 'yu', 'yubar']

        # See https://arxiv.org/abs/2112.10787 for more information about
        # specific couplings and constants

        vev_gev: float = 246
        self.vev: float = vev_gev * GEV

        # Full-SM MS-bar gauge and Yukawa inputs at MZ, Huang and Zhou,
        # Tables II, III and V: https://inspirehep.net/files/71930f37457ebd10742e55c2b36531ea
        # The tabulated masses are in GeV and obey m_f = y_f * v_F / sqrt(2),
        # with v_F = 246 GeV.
        self.yu: np.matrix = np.sqrt(2) * np.diag([0.00123, 0.620, 168.26]) / vev_gev
        self.yd: np.matrix = np.sqrt(2) * np.diag([0.00267, 0.05316, 2.839]) / vev_gev
        self.yl: np.matrix = np.sqrt(2) * np.diag([0.00048307, 0.101766, 1.72856]) / vev_gev

        # Barred matrices are component-wise conjugates, held as separate arrays
        self.yubar: np.matrix = np.conjugate(self.yu)
        self.ydbar: np.matrix = np.conjugate(self.yd)
        self.ylbar: np.matrix = np.conjugate(self.yl)

        self.mZ: float = MZ
        self.mH: float = 125.10 * GEV

        self.g1: float = 0.357254  # U(1)_Y, without GUT normalization
        self.g2: float = 0.65100   # SU(2)_L
        self.g3: float = 1.2104    # SU(3)_c

        # Constants
        self.iCPV: int = 1  # Fixes \eps_{0123} sign convention
        self.onelooporder: int = 1  # Dummy variable to identify 1-loop part
        self.aEV: float = 1  # Evanescent parameter
        self.invepsilonbar: float = 0
        self.epsilonbar: float = 0
        self.xRP: float = 1 # Reading-point parameter

    @property
    def muH(self):
        return np.sqrt(self.mH**2 / 2)

    @property
    def lam(self):
        return self.mH**2 / (2 * self.vev**2)

    @property
    def mu(self):
        return self.scale

    def kronecker_delta(self, i, j):
        if i == j:
            return 1
        return 0

    def get_mass(self):
        return getattr(self, f"M{self.name}")

    def set_mass(self, val):
        setattr(self, f"M{self.name}", val)

    def just_running(self, coeff: str, args: list[int] = []):
        """Return just the non-analytic part of the coefficient `coeff`. The
        method is called with the arguments passed to the list `args`.

        """
        # Save the initial value of the scale
        initial_scale = self.scale

        # Set scale to the mass to kill running contribution
        self.scale = getattr(self, f"M{self.name}")

        # Extract method into a function
        func = getattr(self, coeff)

        # Isolate the matching part
        matching = func(*args)

        # Reset scale to initial value
        self.scale = initial_scale

        # Subtract matching from total result to isolate the running
        running = func(*args) - matching

        return running

    def ratio_tree_loop(self, coeff: str, args: list[int] = []):
        """Returns the ratio between the tree-level part and loop-level part
        `coeff` for flavour indices listed in `args`.

        """
        # Set this here to set it back at the end
        original_onelooporder = self.onelooporder

        # Extract the method into a function
        func = getattr(self, coeff)

        # Get total result
        self.onelooporder = 1
        total_result = func(*args)

        # Get tree-level result
        self.onelooporder = 0
        tree_level_result = func(*args)

        # Isolate just the loop-level result
        loop_level_result = total_result - tree_level_result

        # Set onelooporder back to what it was
        self.onelooporder = original_onelooporder

        return tree_level_result / loop_level_result

    def coeff_dict(self):
        """Return baryon-conserving Warsaw coefficients in LSME TeV units.

        Keys are WCxf Warsaw operator names, but values use the package's own
        units: dimension-six coefficients in TeV**-2 and the dimension-five
        Weinberg coefficients in TeV**-1, as returned by the alpha... methods.
        This is not WCxf input; use to_wcxf for Wilson and other WCxf tools.

        The unsupported baryon-violating operator families are omitted,
        independently of WCxf's sector names. Other missing coefficient
        methods raise AttributeError.
        """
        sectors = wcxf.Basis["SMEFT", "Warsaw"].sectors
        coeff_dict_ = {}
        for ops in sectors.values():
            for op in ops:
                if op.split("_", 1)[0] in BARYON_VIOLATING_WARSAW_OPERATORS:
                    continue
                coeff_dict_[op] = self._get_warsaw_coeff(op)

        return coeff_dict_

    def to_wcxf(self):
        """Return the coefficients of coeff_dict as a WCxf SMEFT Warsaw object.

        Values are converted to WCxf units, GeV**-2 at dimension six and
        GeV**-1 for the Weinberg coefficients, and the scale to GeV. Use
        wilson.Wilson.from_wc(result.to_wcxf()) for running, or .dump() to
        write a WCxf file.
        """
        values = {
            op: self._wcxf_unit_factor(op) * coeff
            for op, coeff in self.coeff_dict().items()
        }
        return wcxf.WC(
            eft="SMEFT",
            basis="Warsaw",
            scale=self.scale / GEV,
            values=wcxf.WC.dict2values(values),
        )

    @staticmethod
    def _wcxf_unit_factor(op: str):
        """Return the TeV-to-GeV factor for one Warsaw coefficient."""
        return GEV if op.split("_", 1)[0] == "llphiphi" else GEV**2

    def _get_warsaw_coeff(self, op: str):
        """Map one Warsaw operator to its coefficient in TeV units."""
        name, *flavour = op.split("_")
        zero_indexed_flavour_indices = [int(i)-1 for i in flavour[0]] if flavour else []
        name = name.replace("phi", "H").replace("tilde", "t")
        wcxf_mme_map = {
            "G": "3G",
            "Gt": "3Gt",
            "W": "3W",
            "Wt": "3Wt",
            "HWtB": "HWBt",
        }

        if name in wcxf_mme_map:
            name = wcxf_mme_map[name]

        func = getattr(self, f"alphaO{name}" if name != "llHH" else "alphaWeinberg")
        return func(*zero_indexed_flavour_indices)

    def exotic_param_dims_dict(self):
        return {
            param_string: np.array(getattr(self, param_string)).shape
            for param_string in self.exotic_params
        }
