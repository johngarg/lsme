#!/usr/bin/env python3

from functools import lru_cache
import numpy as np

from wilson import wcxf

from .sm_running import SM_PARAMETER_NAMES, sm_parameters

GEV = 1e-3
IDENTICAL_CURRENT_OPERATORS = ("ll", "qq1", "qq3", "uu", "dd", "ee")
BARYON_VIOLATING_WARSAW_OPERATORS = frozenset({"duql", "duue", "qqql", "qque"})
# MatchMaker uses D = partial - i g A; WCxf Warsaw uses partial + i g A.
# A -> -A flips every operator containing an odd number of field strengths.
ODD_FIELD_STRENGTH_OPERATORS = frozenset({
    "uB", "uW", "uG", "dB", "dW", "dG", "eB", "eW",
    "G", "Gtilde", "W", "Wtilde",
})

class GenericMatchingResult:
    """Class for storing MatchMaker matching results. Model parameters, direct
    coefficient methods and coeff_dict use TeV units and MatchMaker
    conventions; to_wcxf converts to WCxf conventions and GeV units.

    Scale is initialised to 1 TeV and must lie between MZ and 1000 TeV. The
    SM parameters g1, g2, g3, lam, muH, yu, yd, yl and the barred Yukawas are
    Wilson's pure-SM MS-bar values at the scale, run at one loop from MZ (see
    lsme.numeric.sm_running for the conventions; muH is imaginary). Setting
    scale reloads all of them, overwriting manual values, so set custom SM
    inputs after the scale. Manual SM inputs affect matching only; to_wcxf
    does not pass them to Wilson.

    """

    def __init__(self, name: str, scale: float = 1.0):
        self.name = name
        self.exotic_params = {}
        self.sm_params = list(SM_PARAMETER_NAMES)

        # See https://arxiv.org/abs/2112.10787 for more information about
        # specific couplings and constants
        self.scale = scale

        # Constants
        self.iCPV: int = 1  # Fixes \eps_{0123} sign convention
        self.onelooporder: int = 1  # Dummy variable to identify 1-loop part
        self.aEV: float = 1  # Evanescent parameter
        self.invepsilonbar: float = 0
        self.epsilonbar: float = 0
        self.xRP: float = 1 # Reading-point parameter

    @property
    def scale(self):
        return self._scale

    @scale.setter
    def scale(self, value):
        # Look up first, so an invalid scale leaves the object unchanged.
        parameters = sm_parameters(value)
        self._scale = float(value)
        for name, parameter in parameters.items():
            setattr(self, name, parameter)

    @property
    def mu(self):
        return self._scale

    def kronecker_delta(self, i, j):
        if i == j:
            return 1
        return 0

    def get_mass(self):
        return getattr(self, f"M{self.name}")

    def set_mass(self, val):
        setattr(self, f"M{self.name}", val)

    def just_running(self, coeff: str, args: list[int] = [], reference_scale=None):
        """Return the explicit-log part of the coefficient `coeff` at fixed
        couplings: its value at the current scale minus its value at
        `reference_scale` (default: the model mass), with the SM couplings
        held at their current values. The method is called with the arguments
        passed to the list `args`. Models with several masses must pass
        `reference_scale`; one reference scale cannot remove every mass
        logarithm.

        """
        if reference_scale is None:
            reference_scale = getattr(self, f"M{self.name}")
        reference_scale = float(reference_scale)
        # No table lookup happens here, so only finiteness and positivity
        # are required, not the table's range.
        if not np.isfinite(reference_scale) or reference_scale <= 0:
            raise ValueError(f"Reference scale must be finite and positive, got {reference_scale}")
        func = getattr(self, coeff)
        initial_scale = self._scale
        try:
            # Assign _scale directly so the couplings stay frozen
            self._scale = reference_scale
            matching = func(*args)
        finally:
            self._scale = initial_scale

        return func(*args) - matching

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
        """Return baryon-conserving MatchMaker coefficients with Warsaw keys.

        Keys are WCxf Warsaw operator names, but values retain MatchMaker's
        signs, full-flavour-sum normalisation and the package's own units:
        dimension-six coefficients in TeV**-2 and the dimension-five
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
        """Return a WCxf SMEFT Warsaw object with converted conventions and units.

        Sum distinct current orderings for each identical-current operator,
        including Fierz-related orderings for ee. Operators with an odd number
        of field strengths change sign. The raw alpha... methods and coeff_dict
        retain MatchMaker's conventions; see the README.

        Values are converted to WCxf units, GeV**-2 at dimension six and
        GeV**-1 for the Weinberg coefficients, and the scale to GeV. Use
        wilson.Wilson.from_wc(result.to_wcxf()) for running, or .dump() to
        write a WCxf file.
        """
        plan = _wcxf_export_plan()
        values = {}
        for key, value in self.coeff_dict().items():
            method_name, partners, real, factor = plan[key]
            if partners:
                method = getattr(self, method_name)
                value = value + method(*partners[0])
                if len(partners) == 3:
                    # Add the Fierz pair together, as in the tensor projection.
                    value += method(*partners[1]) + method(*partners[2])
            if real:
                value = value.real
            values[key] = factor * value
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

    @staticmethod
    @lru_cache(maxsize=None)
    def _warsaw_method(op: str):
        """Return the coefficient method name and zero-indexed flavour
        indices for one Warsaw operator."""
        name, *flavour = op.split("_")
        zero_indexed_flavour_indices = tuple(int(i) - 1 for i in flavour[0]) if flavour else ()
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

        method = f"alphaO{name}" if name != "llHH" else "alphaWeinberg"
        return method, zero_indexed_flavour_indices

    def _get_warsaw_coeff(self, op: str):
        """Map one Warsaw operator to its coefficient in TeV units."""
        method, indices = self._warsaw_method(op)
        return getattr(self, method)(*indices)

    def exotic_param_dims_dict(self):
        return {
            param_string: np.array(getattr(self, param_string)).shape
            for param_string in self.exotic_params
        }


@lru_cache(maxsize=1)
def _wcxf_export_plan():
    """Map Warsaw keys to index permutations and conventions shared by all models."""
    plan = {}
    for operators in wcxf.Basis["SMEFT", "Warsaw"].sectors.values():
        for key, metadata in operators.items():
            family = key.partition("_")[0]
            if family in BARYON_VIOLATING_WARSAW_OPERATORS:
                continue
            method, indices = GenericMatchingResult._warsaw_method(key)
            partners = ()
            if family in IDENTICAL_CURRENT_OPERATORS:
                p, r, s, t = indices
                orbit = [indices, (s, t, p, r)]
                if family == "ee":
                    orbit.extend([(s, r, p, t), (p, t, s, r)])
                # coeff_dict already evaluates the first, canonical ordering.
                partners = tuple(dict.fromkeys(orbit))[1:]
            factor = GenericMatchingResult._wcxf_unit_factor(key)
            if family in ODD_FIELD_STRENGTH_OPERATORS:
                factor = -factor
            plan[key] = (method, partners, metadata.get("real", False), factor)
    return plan
