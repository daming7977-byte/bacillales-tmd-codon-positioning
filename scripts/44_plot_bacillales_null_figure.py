from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

NULL_FILE = Path(
    "supplementary/bacillales/"
    "Table_B8_synonymous_permutation_null.tsv"
)

SUMMARY_FILE = Path(
    "supplementary/bacillales/"
    "Table_B9_final_summary.tsv"
)

OUTDIR = Path("figures/bacillales")
OUTDIR.mkdir(parents=True, exist_ok=True)

null_df = pd.read_csv(NULL_FILE, sep="\t")
summary_df = pd.read_csv(SUMMARY_FILE, sep="\t")

anchors = [
    ("start", "TMD start"),
    ("end", "TMD end"),
    ("center", "TMD center"),
]

fig, axes = plt.subplots(
    1, 3,
    figsize=(14, 4.4),
    sharey=True
)

for ax, (anchor, title) in zip(axes, anchors):

    obs = summary_df[
        summary_df["anchor"] == anchor
    ].iloc[0]

    vals = null_df[
        f"{anchor}_median_variance"
    ].dropna()

    ax.hist(vals, bins=30)

    ax.axvline(
        obs["observed_median_variance"],
        linewidth=2.2,
        label="Observed"
    )

    ax.axvline(
        obs["null_median_variance"],
        linewidth=2.2,
        linestyle="--",
        label="Null median"
    )

    ax.set_title(title, fontsize=13)

    ax.set_xlabel(
        "Median cross-species variance (aa²)"
    )

    text = (
        f"Observed = {obs['observed_median_variance']:.2f}\n"
        f"Null median = {obs['null_median_variance']:.2f}\n"
        f"Units = {int(obs['observed_units'])}\n"
        f"P = {obs['empirical_p']:.4f}"
    )

    ax.text(
        0.97,
        0.97,
        text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9.5,
        bbox=dict(
            boxstyle="round,pad=0.35",
            alpha=0.12
        )
    )

axes[0].set_ylabel("Permutation count")

fig.text(
    0.5,
    0.01,
    "Solid vertical line: observed median; "
    "dashed vertical line: permutation-null median",
    ha="center",
    fontsize=9.5
)

fig.tight_layout(
    rect=[0, 0.06, 1, 1]
)

png = OUTDIR / "bacillales_permutation_null_3panel.png"
pdf = OUTDIR / "bacillales_permutation_null_3panel.pdf"

fig.savefig(
    png,
    dpi=300,
    bbox_inches="tight"
)

fig.savefig(
    pdf,
    bbox_inches="tight"
)

plt.close(fig)

print("Saved:")
print(png)
print(pdf)
