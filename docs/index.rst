.. brainpy documentation master file, created by
   sphinx-quickstart on Thu Oct 15 01:53:09 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to brainpy's documentation!
===================================

:mod:`brainpy` is a small Python library implementing the *B* afflingly *R* ecursive
*A* lgorithm for *I* sotopic Patter *N* generation [Dittwald2014]. It includes three implementations,
a pure-Python object oriented implementation, a :title-reference:`Cython` accelerated
version of the object oriented implementation, and a pure :title-reference:`C` implementation,
listed in order of ascending speed. The C implementation is used by default when available.


BRAIN takes an elemental composition represented by any :class:`Mapping`-like Python object
and uses it to compute its aggregated isotopic distribution. All isotopic variants of the same
number of neutrons are collapsed into a single centroid peak, meaning it does not consider
isotopic fine structure.

.. plot::
    :include-source:

    from brainpy import isotopic_variants

    # Generate theoretical isotopic pattern
    peptide = {'H': 53, 'C': 34, 'O': 15, 'N': 7}
    theoretical_isotopic_cluster = isotopic_variants(peptide, npeaks=5, charge=1)
    for peak in theoretical_isotopic_cluster:
        print(peak.mz, peak.intensity)

    # produce a theoretical profile using a gaussian peak shape
    import numpy as np
    mz_grid = np.arange(theoretical_isotopic_cluster[0].mz - 1,
                        theoretical_isotopic_cluster[-1].mz + 1, 0.02)
    intensity = np.zeros_like(grid)
    sigma = 0.002
    for peak in theoretical_isotopic_cluster:
        # Add gaussian peak shape centered around each theoretical peak
        intensity += peak.intensity * np.exp(-(mz_grid - peak.mz) ** 2 / (2 * sigma)
                ) / (np.sqrt(2 * np.pi) * sigma)

    # Normalize profile to 0-100
    intensity = (intensity / intensity.max()) * 100

    # draw the profile
    from matplotlib import pyplot as plt
    plt.plot(grid, intensity)
    plt.xlabel("m/z")
    plt.ylabel("Relative intensity")


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

.. [Dittwald2014]
    Dittwald, P., & Valkenborg, D. (2014). BRAIN 2.0: time and memory complexity improvements in the algorithm for calculating the isotope distribution. Journal of the American Society for Mass Spectrometry, 25(4), 588–94. https://doi.org/10.1007/s13361-013-0796-5
