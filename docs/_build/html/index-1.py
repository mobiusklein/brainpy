from brainpy import isotopic_variants

peptide = {'H': 53, 'C': 34, 'O': 15, 'N': 7}
theoretical_isotopic_cluster = isotopic_variants(peptide, npeaks=5, charge=1)
for peak in theoretical_isotopic_cluster:
    print(peak.mz, peak.intensity)

# produce a theoretical profile using a gaussian peak shape
import numpy as np
grid = np.arange(theoretical_isotopic_cluster[0].mz - 1,
                 theoretical_isotopic_cluster[-1].mz + 1, 0.02)
intensity = np.zeros_like(grid)
sigma = 0.002
for i, mz in enumerate(grid):
    for peak in theoretical_isotopic_cluster:
        intensity[i] += peak.intensity * np.exp(-(mz - peak.mz) ** 2 / (2 * sigma)
                ) / (np.sqrt(2 * np.pi) * sigma)

intensity = (intensity / intensity.max()) * 100

# draw the profile
from matplotlib import pyplot as plt
plt.plot(grid, intensity)
plt.xlabel("m/z")
plt.ylabel("Relative intensity")