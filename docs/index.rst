.. brainpy documentation master file, created by
   sphinx-quickstart on Thu Oct 15 01:53:09 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to brainpy's documentation!
===================================

.. contents:: Table of Contents
    :depth: 3


:mod:`brainpy` is a small Python library implementing the *B* afflingly *R* ecursive
*A* lgorithm for *I* sotopic Patter *N* generation [Dittwald2014]_. It includes three implementations,
a pure-Python object oriented implementation, a :title-reference:`Cython` accelerated
version of the object oriented implementation, and a pure :title-reference:`C` implementation,
listed in order of ascending speed. The C implementation is used by default when available.

BRAIN, implemented in :func:`brainpy.isotopic_variants`, takes an elemental
composition represented by  any :class:`~.collections.abc.Mapping`-like Python object and uses
it to compute its aggregated isotopic distribution. All isotopic variants of the same number
of neutrons are collapsed into a single centroid peak, meaning it does not consider isotopic
fine structure.

.. plot::
    :include-source:

    from brainpy import isotopic_variants

    # Generate theoretical isotopic pattern
    peptide = {'H': 53, 'C': 34, 'O': 15, 'N': 7}
    theoretical_isotopic_cluster = isotopic_variants(peptide, npeaks=5, charge=1)
    for peak in theoretical_isotopic_cluster:
        print(peak.mz, peak.intensity)

    # All following code is to illustrate what brainpy just did.

    # produce a theoretical profile using a gaussian peak shape
    import numpy as np
    mz_grid = np.arange(theoretical_isotopic_cluster[0].mz - 1,
                        theoretical_isotopic_cluster[-1].mz + 1, 0.02)
    intensity = np.zeros_like(mz_grid)
    sigma = 0.002
    for peak in theoretical_isotopic_cluster:
        # Add gaussian peak shape centered around each theoretical peak
        intensity += peak.intensity * np.exp(-(mz_grid - peak.mz) ** 2 / (2 * sigma)
                ) / (np.sqrt(2 * np.pi) * sigma)

    # Normalize profile to 0-100
    intensity = (intensity / intensity.max()) * 100

    # draw the profile
    from matplotlib import pyplot as plt
    plt.plot(mz_grid, intensity)
    plt.xlabel("m/z")
    plt.ylabel("Relative intensity")


Installing
----------
:mod:`brainpy` has three implementations, a pure Python implementation, a Cython translation
of that implementation, and a pure C implementation that releases the :title-reference:`GIL`.

To install from a package index, you will need to have a C compiler appropriate to your Python
version to build these extension modules. Additionally, there are prebuilt wheels for Windows
available on `PyPI <https://pypi.org/project/brain-isotopic-distribution/>`_:

.. code-block:: sh

    $ pip install brain-isotopic-distribution

To build from source, in addition to a C compiler you will also need to install a recent version
of `Cython <https://pypi.org/project/Cython/>`_ to transpile C code.


Usage
-----

:mod:`brainpy` provides a single top-level function for taking a :class:`~collections.abc.Mapping`-like object
defining a elemental composition and generating an isotopic pattern from it, :func:`brainpy.isotopic_variants`.


Module Reference
----------------

.. automodule:: brainpy

    .. autofunction:: isotopic_variants

    .. autofunction:: max_variants

    .. autofunction:: calculate_mass

    .. autofunction:: parse_formula

    .. autoclass:: PyComposition
        :members:

Supporting Objects
******************

    .. autoclass:: Peak

    .. autoclass:: IsotopicDistribution

        .. automethod:: aggregated_isotopic_variants



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Bibliography
============

.. [Dittwald2014]
    Dittwald, P., & Valkenborg, D. (2014). BRAIN 2.0: time and memory complexity improvements in the algorithm for calculating the isotope distribution. Journal of the American Society for Mass Spectrometry, 25(4), 588–94. https://doi.org/10.1007/s13361-013-0796-5
