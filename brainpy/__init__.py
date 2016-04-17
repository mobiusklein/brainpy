'''
A Python Implementation of the Baffling Recursive Algorithm for Isotopic cluster distributioN
'''
import os

from brainpy import (isotopic_variants, IsotopicDistribution, periodic_table,
                     max_variants, calculate_mass, neutral_mass, mass_charge_ratio,
                     PROTON, _has_c, Peak)


def get_include():
    return os.path.join(__path__[0], "_c")

if _has_c:
    from brainpy import _IsotopicDistribution

__author__ = "Joshua Klein & Han Hu"
