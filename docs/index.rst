.. brainpy documentation master file, created by
   sphinx-quickstart on Thu Oct 15 01:53:09 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to brainpy's documentation!
===================================

:mod:`brainpy` is a small Python library implementing the *B*afflingly *R*ecursive
*A*lgorithm for *I*sotopic Patter*N* generation. It includes three implementations,
a pure-Python object oriented implementation, a :title-reference:`Cython` accelerated
version of the object oriented implementation, and a pure :title-reference:`C` implementation,
listed in order of ascending speed. The C implementation is used by default when available.


BRAIN takes an elemental composition represented by any :class:`Mapping`-like Python object
and uses it to compute its aggregated isotopic distribution. All isotopic variants of the same
number of neutrons are collapsed into a single peak, meaning it does not consider isotopic fine
structure.

.. code-block:: python

    from brainpy import isotopic_variants
    
    leucine = dict(C=6, H=13, O=2, N=1)
    theoretical_isotopic_cluster = isotopic_variants(leucine, n_peaks=6, charge=1)
    for peak in theoretical_isotopic_cluster:
        print(peak.mz, peak.intensity)

.. automodule:: brainpy
    
    .. autofunction:: isotopic_variants
    
    .. autofunction:: max_variants
    
    .. autofunction:: calculate_mass
    
    .. autoclass:: Peak

    .. autoclass:: IsotopicDistribution

        .. automethod:: aggregated_isotopic_variants

    .. autofunction:: parse_formula

    .. autoclass:: PyComposition


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

