from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from paths import results_path

# Settings
delta_lower = 20
delta_upper = 50

source_folder = results_path / "tables_first_half" / f"{delta_lower}_to_{delta_upper}_degree_window"
save_folder = results_path / "figures_first_half" / "paired_t_test_res_vector"
save_folder.mkdir(parents=True, exist_ok=True)

# Read and combine the results, ordered by lag
files = sorted(source_folder.glob("paired_bias_stats_N_*_folded.csv"))

if not files:
    raise FileNotFoundError(f"No results found in {source_folder}")

data = pd.concat([pd.read_csv(file, sep=";") for file in files], ignore_index=True)
data = data.loc[data["Lag"].between(2, 9)].sort_values("Lag")

if data["Lag"].duplicated().any():
    raise ValueError("Multiple results were found for the same lag.")

lags = data["Lag"].to_numpy()
mean_difference = data["MeanDifference"].to_numpy()
ci_lower = data["CILower"].to_numpy()
ci_upper = data["CIUpper"].to_numpy()

# Matplotlib expects distances from the mean to each confidence limit
error_bars = np.vstack([mean_difference - ci_lower, ci_upper - mean_difference])

# Plot appearance
plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "font.family": "DejaVu Sans",
    "pdf.fonttype": 42,
})

fig, ax = plt.subplots(figsize=(10.5, 6.4))
fig.subplots_adjust(left=0.12, right=0.97, top=0.80, bottom=0.25)

fig.text(0.12, 0.93, "Bias difference by trial lag", fontsize=20, weight="bold", color="#203447")
fig.text(0.12, 0.875, f"Aligned − non-aligned resultant history  |  {delta_lower}° < |delta| < {delta_upper}°", fontsize=12, color="#516274")

ax.axhline(0, color="#64748b", linestyle="--", linewidth=1.1)

ax.errorbar(
    lags,
    mean_difference,
    yerr=error_bars,
    fmt="o",
    color="#197f88",
    markersize=7,
    markeredgecolor="white",
    markeredgewidth=1,
    elinewidth=2,
    capsize=5,
    zorder=3,
)

ax.set_xticks(lags)
ax.set_xticklabels([f"N−{int(lag)}" for lag in lags])
ax.set_xlabel("Trial lag", labelpad=12)
ax.set_ylabel("Mean folded-error difference (°)", labelpad=12)
ax.set_xlim(lags.min() - 0.6, lags.max() + 0.6)

# Leave additional space above the confidence intervals for p-value labels
plot_range = max(ci_upper.max(), 0) - min(ci_lower.min(), 0)
padding = max(plot_range * 0.15, 0.1)
ax.set_ylim(min(ci_lower.min(), 0) - padding, max(ci_upper.max(), 0) + 1.5 * padding)

ax.set_axisbelow(True)
ax.grid(axis="y", color="#e5e9ed", linewidth=0.8)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_color("#aeb8c2")
ax.spines["bottom"].set_color("#aeb8c2")

# Add the saved p-value above each confidence interval
for row in data.itertuples():
    label = f"p = {row.PValue:.3f}" if row.PValue >= 0.001 else f"p = {row.PValue:.1e}"
    ax.annotate(label, (row.Lag, row.CIUpper), xytext=(0, 12), textcoords="offset points", ha="center", fontsize=9, color="#405164")

# Figure captions
n_min = int(data["NPairs"].min())
n_max = int(data["NPairs"].max())
n_label = str(n_min) if n_min == n_max else f"{n_min}–{n_max}"

fig.text(0.12, 0.135, "Points: mean paired differences. Bars: 95% confidence intervals from the saved paired t-tests.", fontsize=10, color="#516274")
fig.text(0.12, 0.095, f"Positive values indicate greater folded bias in aligned trials. {n_label} observer–condition pairs per lag.", fontsize=10, color="#516274")
fig.text(0.12, 0.055, "Exploratory: repeated-participant dependence is not modelled; p-values are unadjusted across lags.", fontsize=9, color="#516274")

# Save before displaying
figure_name = f"paired_bias_{delta_lower}_to_{delta_upper}"

fig.savefig(save_folder / f"{figure_name}.png", dpi=300, bbox_inches="tight")
fig.savefig(save_folder / f"{figure_name}.pdf", bbox_inches="tight")

plt.show()
plt.close(fig)