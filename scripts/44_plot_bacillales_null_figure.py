from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

NULL_FILE = Path("supplementary/bacillales/Table_B8_synonymous_permutation_null.tsv")
SUMMARY_FILE = Path("supplementary/bacillales/Table_B9_final_summary.tsv")
OUTDIR = Path("figures/bacillales")
OUTDIR.mkdir(parents=True, exist_ok=True)

null_df = pd.read_csv(NULL_FILE, sep="\t")
summary_df = pd.read_csv(SUMMARY_FILE, sep="\t")

anchors = ["start", "end", "center"]

fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))

for ax, anchor in zip(axes, anchors):
    obs = summary_df.loc[summary_df["anchor"] == anchor].iloc[0]

    null_col = f"{anchor}_median_variance"
    vals = null_df[null_col].dropna()

    ax.hist(vals, bins=30)
    ax.axvline(obs["observed_median_variance"], linewidth=2, label="Observed")
    ax.axvline(obs["null_median_variance"], linewidth=2, linestyle="--", label="Null median")

    ax.set_title(anchor.capitalize())
    ax.set_xlabel("Median variance")
    ax.set_ylabel("Permutation count")

    text = (
        f"Observed = {obs['observed_median_variance']:.2f}\n"
        f"Null median = {obs['null_median_variance']:.2f}\n"
        f"Units = {int(obs['observed_units'])}\n"
        f"P = {obs['empirical_p']:.4f}"
    )
    ax.text(
        0.98, 0.98, text,
        transform=ax.transAxes,
        ha="right", va="top",
        bbox=dict(boxstyle="round", alpha=0.15)
    )

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False)
fig.suptitle("Independent Bacillales analysis: synonymous-permutation null distributions", y=1.03)
fig.tight_layout()

png = OUTDIR / "bacillales_permutation_null_3panel.png"
pdf = OUTDIR / "bacillales_permutation_null_3panel.pdf"

fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")

print("Saved:")
print(png)
print(pdf)
