"""
Build a figure showing the BayesTraits posterior distributions of T6SS transition rates in motile and non-motile lineages.

Half-page target size: ~130 mm × 90 mm
Font: Helvetica throughout

Each panel contains two rows:
  A  T6SS gain rate. Top: q34 (motile). Bottom: q12 (non-motile).
  B  T6SS loss rate. Top: q43 (motile). Bottom: q21 (non-motile).

Input file required in the same folder:
  bayestraits_rates.tsv

Expected columns:
  q12, q21, q34, q43
"""

import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Font check: use Helvetica only
# ---------------------------------------------------------------------------
available_fonts = sorted(set(f.name for f in font_manager.fontManager.ttflist))
if "Helvetica" not in available_fonts:
    print(
        "WARNING: Helvetica was not found on this system. "
        "Matplotlib may substitute another font. For final submission, "
        "open the PDF in Illustrator and set all text to Helvetica."
    )

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"],
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "xtick.major.width": 0.9,
    "ytick.major.width": 0.9,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "legend.frameon": False,
    "legend.fontsize": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# Palette
COL_MOT = "#7E70AE"      # motile
COL_NONMOT = "#5B7FAE"   # non-motile

# ---------------------------------------------------------------------------
# Load posterior samples
# ---------------------------------------------------------------------------
input_file = os.path.join(HERE, "bayestraits_rates.tsv")

data = np.genfromtxt(
    input_file,
    delimiter="\t",
    names=True,
)

q12 = data["q12"]   # non-motile T6SS gain
q21 = data["q21"]   # non-motile T6SS loss
q34 = data["q34"]   # motile T6SS gain
q43 = data["q43"]   # motile T6SS loss

# ---------------------------------------------------------------------------
# Helper function to draw one two-row panel
# ---------------------------------------------------------------------------
def draw_two_row_panel(ax, vals_mot, vals_non, xlabel, title,
                       zero_line=False, rng_seed=0):
    rng = np.random.default_rng(rng_seed)

    if zero_line:
        ax.axvline(0, color="#888888", linewidth=0.9, linestyle="--", zorder=1)

    groups = [
        ("motile", vals_mot, COL_MOT, 1.0),
        ("non-motile", vals_non, COL_NONMOT, 0.0),
    ]

    for _, vals, col, y in groups:
        # Violin
        parts = ax.violinplot(
            [vals],
            positions=[y],
            vert=False,
            widths=0.7,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )
        for body in parts["bodies"]:
            body.set_facecolor(col)
            body.set_alpha(0.30)
            body.set_edgecolor(col)
            body.set_linewidth(0.8)

        # Jittered strip
        y_jit = y + rng.uniform(-0.13, 0.13, size=len(vals))
        ax.scatter(
            vals,
            y_jit,
            s=4,
            color=col,
            alpha=0.30,
            edgecolor="none",
            zorder=2,
        )

        # IQR bar and median dot
        q25, q50, q75 = np.percentile(vals, [25, 50, 75])
        ax.plot([q25, q75], [y, y], color=col, linewidth=2.6, zorder=4)
        ax.scatter(
            [q50],
            [y],
            color="white",
            edgecolor=col,
            s=42,
            linewidth=1.5,
            zorder=5,
        )

    # Median annotations above each violin
    for _, vals, col, y in groups:
        med = np.median(vals)
        ax.text(
            med,
            y + 0.46,
            f"median = {med:+.2f}",
            fontsize=8,
            color=col,
            ha="center",
            va="center",
            bbox=dict(
                boxstyle="round,pad=0.25",
                facecolor="white",
                edgecolor="none",
            ),
        )

    # Axes
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["non-motile", "motile"])
    ax.set_ylim(-0.5, 1.7)
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=10, pad=6)

    # X-axis padding
    all_vals = np.concatenate([vals_mot, vals_non])
    xmin, xmax = all_vals.min(), all_vals.max()
    span = xmax - xmin
    pad = 0.10 * span if span > 0 else 1

    if zero_line:
        ax.set_xlim(min(xmin, 0) - pad, max(xmax, 0) + pad)
    else:
        ax.set_xlim(xmin - pad, xmax + pad)


def panel_label(ax, label, x=-0.06, y=1.02):
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        va="bottom",
        ha="left",
    )

# ---------------------------------------------------------------------------
# Build half-page figure
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(
    1,
    2,
    figsize=(5.1, 3.5),  # approximately 130 mm × 89 mm
    gridspec_kw=dict(
        wspace=0.35,
        left=0.10,
        right=0.98,
        top=0.85,
        bottom=0.22,
    ),
)

draw_two_row_panel(
    axes[0],
    q34,
    q12,
    xlabel=r"T6SS gain rate ($q_{34}$, $q_{12}$)",
    title="T6SS gain rate",
    rng_seed=1,
)
panel_label(axes[0], "A")

draw_two_row_panel(
    axes[1],
    q43,
    q21,
    xlabel=r"T6SS loss rate ($q_{43}$, $q_{21}$)",
    title="T6SS loss rate",
    rng_seed=2,
)
panel_label(axes[1], "B")

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
out_pdf = os.path.join(HERE, "bayestraits_posterior_halfpage_helvetica.pdf")
out_png = os.path.join(HERE, "bayestraits_posterior_halfpage_helvetica.png")

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, dpi=600, bbox_inches="tight")

print(f"wrote {out_pdf}")
print(f"wrote {out_png}")

# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------
def summarise(name, vals):
    q025, q25, q50, q75, q975 = np.percentile(vals, [2.5, 25, 50, 75, 97.5])
    print(
        f"{name}: median={q50:+.3f}, "
        f"50% CI [{q25:+.3f}, {q75:+.3f}], "
        f"95% CI [{q025:+.3f}, {q975:+.3f}]"
    )

print()
summarise("q34 (motile gain)    ", q34)
summarise("q12 (non-motile gain)", q12)
summarise("q43 (motile loss)    ", q43)
summarise("q21 (non-motile loss)", q21)
print()
print(f"P(q34 > q12) = {(q34 > q12).mean():.4f}")
print(f"P(q43 > q21) = {(q43 > q21).mean():.4f}")
