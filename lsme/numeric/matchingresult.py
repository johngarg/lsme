#!/usr/bin/env python3

from dataclasses import dataclass
from typing import List
import numpy as np

from wilson import wcxf

GEV = 1e-3
MZ = 91.1876 * GEV

class GenericMatchingResult:
    """Class for storing MatchMaker matching results. All dimensionful quantities
    are in TeV. Scale is initialised to 1 TeV.

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

        # SM couplings and parameters
        # Values taken from https://inspirehep.net/files/71930f37457ebd10742e55c2b36531ea
        self.yu: np.matrix = np.diag([0.00123, 0.620, 168.26]) / vev_gev
        self.yd: np.matrix = np.diag([0.00267, 0.05316, 2.839]) / vev_gev
        self.yl: np.matrix = np.diag([0.00048307, 0.101766, 1.72856]) / vev_gev

        self.yubar: np.matrix = self.yu
        self.ydbar: np.matrix = self.yd
        self.ylbar: np.matrix = self.yl

        self.mZ: float = MZ
        self.mH: float = 125.10 * GEV

        self.g1: float = 0.65100
        self.g2: float = 0.357254
        self.g3: float = 1.204

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
        sectors = wcxf.Basis["SMEFT", "Warsaw"].sectors
        coeff_dict_ = {}
        for sector, ops in sectors.items():
            if sector == "dB=dL=1":
                print(
                    f"Skipping operators in {sector} sector, since these are not implemented yet."
                )
                continue

            for op in ops.keys():
                coeff_dict_[op] = self._get_wc_in_wcxf(op)

        return coeff_dict_

    def _get_wc_in_wcxf(self, op: str):
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
