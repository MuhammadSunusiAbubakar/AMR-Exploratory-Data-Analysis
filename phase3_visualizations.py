import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Create output directory for high-DPI figures
os.makedirs("figures", exist_ok=True)

# Load cleaned dataset
df = pd.read_csv("cleaned_amr_data.csv")

# Global decluttered publication styling (Unit 5: Data-Ink Ratio)
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial"],
    "axes.edgecolor": "#CCCCCC",
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

ast_cols = [
    "AMX/AMP", "AMC", "CZ", "FOX", "CTX/CRO", "IPM", "GEN", "AN",
    "Acide nalidixique", "ofx", "CIP", "C", "Co-trimoxazole", "Furanes", "colistine"
]

# ==============================================================================
# FIGURE 1: Pathogen Prevalence (Horizontal Bar Chart with Pre-Attentive Accent)
# ==============================================================================
species_counts = df["Species"].value_counts()
total_known = species_counts.sum()
ecoli_pct = (species_counts["Escherichia coli"] / total_known) * 100

fig, ax = plt.subplots(figsize=(9, 5.2))
colors = ["#1F4E78" if sp == "Escherichia coli" else "#B8C4CC" for sp in species_counts.index]

bars = ax.barh(species_counts.index[::-1], species_counts.values[::-1], color=colors[::-1], height=0.68)
ax.set_xlabel("Number of Clinical Isolates (n)", fontsize=10, fontweight="bold", color="#333333")
ax.set_title(f"Escherichia coli Dominates Clinical Isolates ({ecoli_pct:.1f}% of Identified Samples)",
             fontsize=13, fontweight="bold", loc="left", pad=18, color="#111111")
ax.text(0, 1.02, f"Distribution of standardized bacterial species across n = {total_known:,} identified isolates",
        transform=ax.transAxes, fontsize=9.5, color="#555555")

for bar, count in zip(bars, species_counts.values[::-1]):
    pct = (count / total_known) * 100
    ax.text(count + 80, bar.get_y() + bar.get_height()/2,
            f"{count:,} ({pct:.1f}%)", va="center", ha="left", fontsize=9, color="#222222")

ax.set_xlim(0, species_counts.max() * 1.22)
ax.xaxis.grid(True, linestyle="--", alpha=0.4, color="#CCCCCC")
ax.set_axisbelow(True)
plt.tight_layout()
plt.savefig("figures/fig1_species_prevalence.png", dpi=300, bbox_inches="tight")
plt.close()

# ==============================================================================
# FIGURE 2: Distribution of Isolate Resistance Rate (Histogram + Smoothed KDE)
# ==============================================================================
fig, ax = plt.subplots(figsize=(8.5, 4.8))
data_rate = df["Resistance_Rate_Pct"].dropna()
mean_r = data_rate.mean()
med_r = data_rate.median()

sns.histplot(data_rate, bins=14, kde=True, kde_kws={"bw_adjust": 2.0},
             color="#4A89B8", edgecolor="white", alpha=0.65,
             line_kws={"linewidth": 2.2, "color": "#1F4E78"}, ax=ax)

ax.axvline(mean_r, color="#D95F02", linestyle="--", linewidth=1.8)
ax.axvline(med_r, color="#1B9E77", linestyle="-", linewidth=1.8)

ymax = ax.get_ylim()[1]
label_box = dict(facecolor="white", alpha=0.88, edgecolor="none", pad=2)
ax.text(mean_r - 1.5, ymax * 0.92, f"Mean: {mean_r:.1f}%", color="#D95F02",
        ha="right", fontweight="bold", fontsize=9.5, bbox=label_box)
ax.text(med_r + 1.5, ymax * 0.92, f"Median: {med_r:.1f}%", color="#1B9E77",
        ha="left", fontweight="bold", fontsize=9.5, bbox=label_box)

ax.set_xlabel("Antibiotic Resistance Rate per Isolate (%)", fontsize=10, fontweight="bold", color="#333333")
ax.set_ylabel("Number of Isolates", fontsize=10, fontweight="bold", color="#333333")
ax.set_title("Isolate Resistance Rates Cluster Around 32% With Mild Right Skew",
             fontsize=12.5, fontweight="bold", loc="left", pad=18, color="#111111")
ax.text(0, 1.02, f"Histogram (14 bins) and smoothed KDE across n = {len(data_rate):,} tested isolates",
        transform=ax.transAxes, fontsize=9.5, color="#555555")
ax.yaxis.grid(True, linestyle="--", alpha=0.4, color="#CCCCCC")
ax.set_axisbelow(True)
plt.tight_layout()
plt.savefig("figures/fig2_resistance_distribution.png", dpi=300, bbox_inches="tight")
plt.close()

# ==============================================================================
# FIGURE 3: Resistance Rate by Bacterial Species (Sorted Box + Violin Plot)
# ==============================================================================
df_sp = df.dropna(subset=["Species", "Resistance_Rate_Pct"])
order_sp = df_sp.groupby("Species")["Resistance_Rate_Pct"].median().sort_values(ascending=False).index

fig, ax = plt.subplots(figsize=(9.5, 5.5))
sns.violinplot(data=df_sp, y="Species", x="Resistance_Rate_Pct", order=order_sp,
               inner=None, cut=0, color="#D9E6F2", linewidth=0.8, ax=ax)
sns.boxplot(data=df_sp, y="Species", x="Resistance_Rate_Pct", order=order_sp,
            width=0.25, color="#1F4E78", fliersize=2.5,
            boxprops={"alpha": 0.85}, medianprops={"color": "#FFDD44", "linewidth": 2}, ax=ax)

ax.set_xlim(-2, 95)
ax.set_xlabel("Antibiotic Resistance Rate (%)", fontsize=10, fontweight="bold", color="#333333")
ax.set_ylabel("")
ax.set_title("Multi-Drug Resistance Burden Varies Substantially Across Pathogens",
             fontsize=12.5, fontweight="bold", loc="left", pad=18, color="#111111")
ax.text(0, 1.02, "Combined violin (density, clipped at 0%) and box plots (IQR & median in yellow), sorted by median",
        transform=ax.transAxes, fontsize=9.5, color="#555555")
ax.xaxis.grid(True, linestyle="--", alpha=0.4, color="#CCCCCC")
ax.set_axisbelow(True)
plt.tight_layout()
plt.savefig("figures/fig3_species_resistance_violin_box.png", dpi=300, bbox_inches="tight")
plt.close()

# ==============================================================================
# FIGURE 4: Drug-Level Resistance Profile Across All 15 Antibiotics
# ==============================================================================
drug_res = {}
for col in ast_cols:
    valid = df[col].dropna()
    drug_res[col] = (valid == "R").mean() * 100

drug_series = pd.Series(drug_res).sort_values(ascending=True)

fig, ax = plt.subplots(figsize=(9, 6))
bar_colors = ["#C0392B" if val >= 50 else "#6C8EAD" for val in drug_series.values]

bars = ax.barh(drug_series.index, drug_series.values, color=bar_colors, height=0.68)
ax.axvline(50, color="#C0392B", linestyle=":", linewidth=1.2, alpha=0.7)
ax.text(50.8, 0.5, "50% Critical\nResistance Threshold", color="#C0392B", fontsize=8.5, fontweight="bold")

for bar, val in zip(bars, drug_series.values):
    ax.text(val + 1.0, bar.get_y() + bar.get_height()/2,
            f"{val:.1f}%", va="center", ha="left", fontsize=9, color="#222222")

ax.set_xlim(0, max(drug_series.values) * 1.18)
ax.set_xlabel("Percentage of Resistant Isolates (%)", fontsize=10, fontweight="bold", color="#333333")
ax.set_title("Beta-Lactams Show Critical Resistance (>57%), Contrasting With Non-Beta-Lactams (<20%)",
             fontsize=11.5, fontweight="bold", loc="left", pad=18, color="#111111")
ax.text(0, 1.02, "Percentage of tested isolates classified as Resistant (R) per antibiotic (red >= 50%)",
        transform=ax.transAxes, fontsize=9.5, color="#555555")
ax.xaxis.grid(True, linestyle="--", alpha=0.4, color="#CCCCCC")
ax.set_axisbelow(True)
plt.tight_layout()
plt.savefig("figures/fig4_antibiotic_resistance_rates.png", dpi=300, bbox_inches="tight")
plt.close()

# ==============================================================================
# FIGURE 5: Antibiotic Co-Resistance Correlation Heatmap (Diverging Palette)
# ==============================================================================
df_bin = df[ast_cols].apply(lambda col: col.map({"R": 1, "S": 0, "I": 0}))
corr_matrix = df_bin.corr()

mask = np.triu(np.ones_like(corr_matrix, dtype=bool))

fig, ax = plt.subplots(figsize=(10.5, 8.5))
sns.heatmap(corr_matrix, mask=mask, annot=True, fmt=".2f", cmap="RdBu_r",
            center=0, vmin=-0.10, vmax=0.35, square=True, linewidths=0.5,
            cbar_kws={"shrink": 0.75, "label": "Pearson Correlation (r)"},
            annot_kws={"size": 8}, ax=ax)

ax.set_title("Co-Resistance Correlation Matrix Reveals a Distinct Beta-Lactam Cluster",
             fontsize=13, fontweight="bold", loc="left", pad=20, color="#111111")
ax.text(0, 1.015, "Lower-triangle Pearson correlation (R=1, S/I=0); diverging RdBu palette centered at 0 (scale: -0.10 to 0.35)",
        transform=ax.transAxes, fontsize=9.5, color="#555555")
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.yticks(rotation=0, fontsize=9)
plt.tight_layout()
plt.savefig("figures/fig5_coresistance_heatmap.png", dpi=300, bbox_inches="tight")
plt.close()

print("All 5 polished figures successfully saved to 'figures/'!")