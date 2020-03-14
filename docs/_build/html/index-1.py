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