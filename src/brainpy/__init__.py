'''
A Python Implementation of the Baffling Recursive Algorithm for Isotopic cluster distributioN
'''
import os

from .brainpy import (isotopic_variants, IsotopicDistribution, periodic_table,
                      max_variants, calculate_mass, neutral_mass, mass_charge_ratio,
                      PROTON, _has_c, Peak)

from .composition import parse_formula, PyComposition


SimpleComposition = PyComposition


def get_include():
    """Retrieve the path to compiled C extensions' source files to make linking simple.

    This module contains two variants of the algorithm reimplimented using C and the Python-C API.

    The `_speedup` module is a direct translation of the pure Python implementation using Cython,
    using static typing and `cdef class` versions of the existing classes. As this implementation still
    spends a substantial amount of time in Python-space, it is slower than the option below, but is more
    straight-forward to manipulate from Python.

    The `_c` module is a complete rewrite of the algorithm directly in C, using Python only to load
    mass configuration information and is fully usable. It exports an entrypoint function to Python
    which replaces the :func:`isotopic_variants` function when available. Because almost all of the
    action happens in C here, it's not possible to run individual parts of the process directly from
    Python.
    """
    return os.path.join(__path__[0], "_c")


if _has_c:
    from .brainpy import _IsotopicDistribution

__author__ = "Joshua Klein & Han Hu"


__all__ = [
    "isotopic_variants", "IsotopicDistribution", "periodic_table",
    "max_variants", "calculate_mass", "neutral_mass", "mass_charge_ratio",
    "PROTON", "_has_c", "Peak",

    "parse_formula", "PyComposition", "SimpleComposition",

    "_IsotopicDistribution",

    "get_include"
]
