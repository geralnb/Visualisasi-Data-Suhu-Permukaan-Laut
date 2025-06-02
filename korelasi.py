import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

# Data dari hasil observasi/gambar
data = {
    "Tahun": list(range(2014, 2025)),
    "NDVI_mean": [
        0.78529946, 0.7538026, 0.79785436, 0.7868034, 0.77328485,
        0.79279333, 0.78129935, 0.7912825, 0.86598295, 0.7703251, 0.7703251
    ],
    "NDDI_mean": [
        0.37399829, 0.45232222, 0.34203075, 0.36540276, 0.41209644,
        0.38334647, 0.41377044, 0.39790487, 0.4324071, 0.39188457, 0.39880485
    ],
    "SPL_JJA": [
        27.56666692, 26.27333323, 29.12333361, 26.44666672, 25.77666664,
        24.6899999, 27.46333377, 28.06999969, 28.13333321, 26.12999988, 26.90999985
    ]
}

df = pd.DataFrame(data)

# Hitung korelasi Pearson
ndvi_corr, ndvi_p = pearsonr(df["NDVI_mean"], df["SPL_JJA"])
nddi_corr, nddi_p = pearsonr(df["NDDI_mean"], df["SPL_JJA"])

print(f"Korelasi Pearson NDVI vs SPL: {ndvi_corr:.3f} (p={ndvi_p:.3f})")
print(f"Korelasi Pearson NDDI vs SPL: {nddi_corr:.3f} (p={nddi_p:.3f})")

# Plot
plt.figure(figsize=(12, 5))
plt.suptitle("Korelasi NDVI & NDDI terhadap Suhu JJA (2014–2024)", fontsize=14, fontweight='bold')

# NDVI plot
plt.subplot(1, 2, 1)
sns.regplot(x="NDVI_mean", y="SPL_JJA", data=df, color="green", marker='o')
plt.title("NDVI vs Suhu JJA")
plt.xlabel("NDVI Mean")
plt.ylabel("Suhu Rata-rata JJA (°C)")

# NDDI plot
plt.subplot(1, 2, 2)
sns.regplot(x="NDDI_mean", y="SPL_JJA", data=df, color="red", marker='o')
plt.title("NDDI vs Suhu JJA")
plt.xlabel("NDDI Mean")
plt.ylabel("Suhu Rata-rata JJA (°C)")

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()
