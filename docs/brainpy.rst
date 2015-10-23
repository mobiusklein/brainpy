brainpy package
===============

.. code-block:: python

    from brainpy import isotopic_variants
    
    leucine = dict(C=6, H=13, O=2, N=1)
    theoretical_isotopic_cluster = isotopic_variants(leucine, n_peaks=6, charge=1)
    for peak in theoretical_isotopic_cluster:
        print(peak.mz, peak.intensity)

brainpy.brainpy module
----------------------

.. automodule:: brainpy.brainpy
    
    .. autofunction:: isotopic_variants
    .. autofunction:: max_variants
    .. autofunction:: calculate_mass
    .. autoclass:: Peak

Module contents
---------------

.. automodule:: brainpy
    :members:
    :undoc-members:
    :show-inheritance:    
