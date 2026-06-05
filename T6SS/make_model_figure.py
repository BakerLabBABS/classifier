"""
Build a multi-panel figure for the motility/T6SS
kin-competition model that accompanies the comparative-genomics paper.

Panels:
  A  Conceptual schematic: motility extends the producer's exploitation reach
     L(D) but contracts the relatedness scale xi(D).
  B  Contact-weapon fitness W_C vs effective motility D — robust monotonic
     increase (the clean T6SS prediction).
  C  Diffusible-weapon fitness W_D vs D for several conspecific densities —
     context-dependent (some curves rise, some fall, some peak).
  D  Direct vs kin-mediated components of diffusible-weapon payoff vs D —
     reveals the underlying trade-off.
  E  Heatmap of W_D(motile) - W_D(non-motile) over (toxin range ell,
     conspecific density lambda_s). The sign flips.
  F  Predicted macroevolutionary contrast: distribution of motility-induced
     fitness gain across parameter space — narrow positive for contact
     weapons, broad and zero-straddling for diffusible weapons.
"""

import os
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
import matplotlib.patheffects as pe

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from kin_competition_motility_weapon_model import (
    params,
    exploitation_length,
    relatedness_length,
    relatedness_at_distance,
    conspecific_competition_weight,
    focal_access_weight,
    contact_per_kill_benefit,
    contact_encounter_rate,
    diffusible_kill_distance_pdf,
    diffusible_total_kill_opportunity,
    diffusible_mean_per_kill_benefit,
    contact_fitness,
    diffusible_fitness,
)


# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9.5,
    "axes.labelsize": 10,
    "axes.titlesize": 10.5,
    "axes.titleweight": "regular",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.9,
    "xtick.major.width": 0.9,
    "ytick.major.width": 0.9,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "legend.frameon": False,
    "legend.fontsize": 8.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# Palette
COL_NONMOT = "#2E5E8C"   # deep blue
COL_MOT    = "#D7642C"   # warm orange
COL_DIRECT = "#1F7A4C"   # green
COL_KIN    = "#7B3294"   # purple
COL_T6SS   = "#C42E2E"   # crimson
COL_BACT   = "#5B7FAE"   # slate blue
COL_FOCAL  = "#E69F00"   # gold
COL_KINCELL = "#F2C679"  # light gold
COL_COMP   = "#7E7E7E"   # gray (competitor)
COL_NEUTRAL = "#444444"


D_nm = params["D_nonmotile"]
D_m  = params["D_motile"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def panel_label(ax, label, x=-0.18, y=1.06):
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="bottom", ha="left")


def diffusible_components(D, lambda_s, ell):
    """Return (direct, kin, total_W) components of the diffusible payoff."""
    L = exploitation_length(D)
    rmax = max(10 * ell, 5 * L, 3.0)
    r = np.linspace(1e-6, rmax, 4000)
    pdf = diffusible_kill_distance_pdf(r, ell)
    w_focal = focal_access_weight(r, D)
    C_s = conspecific_competition_weight(D, lambda_s)
    R = relatedness_at_distance(r, D)
    B = params["background_resource_loss_weight"]
    denom = w_focal + C_s + B
    direct_share = w_focal / denom
    kin_share = R * C_s / denom
    mean_direct = np.trapezoid(pdf * direct_share, r)
    mean_kin = np.trapezoid(pdf * kin_share, r)
    K = diffusible_total_kill_opportunity(ell)
    b_D = params["diffusible_benefit_scale"]
    c_D = params["diffusible_cost"]
    return (b_D * K * mean_direct,
            b_D * K * mean_kin,
            b_D * K * (mean_direct + mean_kin) - c_D)


# ---------------------------------------------------------------------------
# Build figure
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(13.5, 8.6))
gs = GridSpec(
    2, 3,
    figure=fig,
    left=0.055, right=0.985, top=0.945, bottom=0.075,
    wspace=0.34, hspace=0.42,
)

ax_A = fig.add_subplot(gs[0, 0])
ax_B = fig.add_subplot(gs[0, 1])
ax_C = fig.add_subplot(gs[0, 2])
ax_D = fig.add_subplot(gs[1, 0])
ax_E = fig.add_subplot(gs[1, 1])
ax_F = fig.add_subplot(gs[1, 2])


# ---------------------------------------------------------------------------
# Panel A — Conceptual schematic
# ---------------------------------------------------------------------------
ax = ax_A
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_aspect("equal")
ax.set_xticks([]); ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)

# Two small "arenas" representing non-motile vs motile worlds
def draw_arena(cx, cy, rad, motile):
    # Arena outline
    ax.add_patch(Circle((cx, cy), rad, facecolor="#FAFAFA",
                        edgecolor="#888888", linewidth=0.9, zorder=1))

    rng = np.random.default_rng(2 if motile else 1)
    # Exploitation length L(D) — for the schematic, use a log-scaled
    # visualization so the non-motile L ring is large enough to read.
    # Numeric values are still shown beneath each arena.
    D = D_m if motile else D_nm
    L = exploitation_length(D)
    xi = relatedness_length(D)
    L_min, L_max = 0.06, 2.0       # spans the values used in the figure
    xi_min, xi_max = 0.24, 1.2
    # Log-scale map onto [0.28, 0.97] * rad — keeps non-motile L visible
    def log_map(x, lo, hi, vis_lo=0.30, vis_hi=0.97):
        t = (np.log(x) - np.log(lo)) / (np.log(hi) - np.log(lo))
        return (vis_lo + (vis_hi - vis_lo) * np.clip(t, 0, 1)) * rad
    L_vis  = log_map(L,  L_min,  L_max)
    xi_vis = log_map(xi, xi_min, xi_max)

    # Relatedness halo (broad gradient) — solid disc with low alpha
    ax.add_patch(Circle((cx, cy), xi_vis, facecolor=COL_FOCAL,
                        alpha=0.16, edgecolor="none", zorder=2))
    # Exploitation reach (dashed)
    ax.add_patch(Circle((cx, cy), L_vis, facecolor="none",
                        edgecolor=COL_FOCAL, linewidth=1.1, linestyle=(0, (3, 2)),
                        zorder=3))

    # Focal cell
    ax.add_patch(Circle((cx, cy), 0.018, facecolor=COL_FOCAL,
                        edgecolor="black", linewidth=0.6, zorder=6))

    # Kin cells — clustered for non-motile, dispersed for motile
    if not motile:
        # Cluster kin near focal
        offsets = rng.normal(0, 0.045, size=(7, 2))
        for ox, oy in offsets:
            ax.add_patch(Circle((cx + ox, cy + oy), 0.013,
                                facecolor=COL_KINCELL, edgecolor="black",
                                linewidth=0.4, zorder=5))
        # Competitors farther away (in distinct patches)
        for _ in range(8):
            theta = rng.uniform(0, 2 * np.pi)
            rr = rng.uniform(rad * 0.55, rad * 0.95)
            x, y = cx + rr * np.cos(theta), cy + rr * np.sin(theta)
            ax.add_patch(Circle((x, y), 0.013, facecolor=COL_COMP,
                                edgecolor="black", linewidth=0.4, zorder=4))
    else:
        # Kin and competitors are mixed
        for _ in range(7):
            theta = rng.uniform(0, 2 * np.pi)
            rr = rng.uniform(0.03, rad * 0.95)
            x, y = cx + rr * np.cos(theta), cy + rr * np.sin(theta)
            ax.add_patch(Circle((x, y), 0.013, facecolor=COL_KINCELL,
                                edgecolor="black", linewidth=0.4, zorder=5))
        for _ in range(8):
            theta = rng.uniform(0, 2 * np.pi)
            rr = rng.uniform(0.03, rad * 0.95)
            x, y = cx + rr * np.cos(theta), cy + rr * np.sin(theta)
            ax.add_patch(Circle((x, y), 0.013, facecolor=COL_COMP,
                                edgecolor="black", linewidth=0.4, zorder=4))
        # Motion lines on focal
        for theta in np.linspace(0, 2 * np.pi, 6, endpoint=False):
            x0 = cx + 0.022 * np.cos(theta)
            y0 = cy + 0.022 * np.sin(theta)
            x1 = cx + 0.04 * np.cos(theta)
            y1 = cy + 0.04 * np.sin(theta)
            ax.plot([x0, x1], [y0, y1], color=COL_FOCAL,
                    linewidth=0.9, solid_capstyle="round", zorder=7)


arena_rad = 0.18
draw_arena(0.26, 0.45, arena_rad, motile=False)
draw_arena(0.74, 0.45, arena_rad, motile=True)

# Arena titles
ax.text(0.26, 0.71, "Non-motile  (low $D$)", ha="center", va="bottom",
        fontsize=9.5, fontweight="bold", color=COL_NONMOT)
ax.text(0.74, 0.71, "Motile  (high $D$)", ha="center", va="bottom",
        fontsize=9.5, fontweight="bold", color=COL_MOT)

# Annotations of L and xi for each arena
ax.text(0.26, 0.205,
        f"$L = {exploitation_length(D_nm):.2f}$    $\\xi = {relatedness_length(D_nm):.2f}$",
        ha="center", va="top", fontsize=8.8, color="#333")
ax.text(0.74, 0.205,
        f"$L = {exploitation_length(D_m):.2f}$    $\\xi = {relatedness_length(D_m):.2f}$",
        ha="center", va="top", fontsize=8.8, color="#333")

# Legend chips at bottom (cell types)
chip_y = 0.085
chips = [
    (0.05, COL_FOCAL,   "focal"),
    (0.17, COL_KINCELL, "kin"),
    (0.27, COL_COMP,    "competitor"),
]
for x, col, lab in chips:
    ax.add_patch(Circle((x, chip_y), 0.014, facecolor=col,
                        edgecolor="black", linewidth=0.5))
    ax.text(x + 0.027, chip_y, lab, va="center", ha="left", fontsize=8.2)

# Right legend: L (dashed) and xi (filled halo)
ax.plot([0.50, 0.575], [chip_y, chip_y], color=COL_FOCAL,
        linestyle=(0, (3, 2)), linewidth=1.4)
ax.text(0.585, chip_y, "exploit  $L(D)$", va="center", fontsize=8.2)
ax.add_patch(Rectangle((0.79, chip_y - 0.018), 0.05, 0.036,
                       facecolor=COL_FOCAL, alpha=0.22, edgecolor="none"))
ax.text(0.855, chip_y, "kin  $\\xi(D)$", va="center", fontsize=8.2)

# Top: title and mechanism descriptor
ax.text(0.5, 0.97,
        "Motility extends $L(D)$ but shrinks $\\xi(D)$",
        ha="center", va="top", fontsize=10.0, color=COL_NEUTRAL)

panel_label(ax, "A", x=-0.04, y=1.02)


# ---------------------------------------------------------------------------
# Panel B — Contact weapon W_C vs D for several lambda_c (competitor density)
# ---------------------------------------------------------------------------
ax = ax_B
D_MAX_PLOT = 0.7   # x-axis upper bound for B, C, D — focuses on the interesting range
D_grid = np.linspace(0.001, D_MAX_PLOT, 240)
lambda_s_default = 0.5
lambda_c_values = [0.3, 1.0, 2.5]
red_cmap = plt.get_cmap("Reds")
red_shades = [red_cmap(0.45), red_cmap(0.65), red_cmap(0.88)]

ax.axhline(0, color="#BBBBBB", linewidth=0.7, zorder=1)
for lc, col in zip(lambda_c_values, red_shades):
    W_C_curve = np.array([contact_fitness(D, lambda_s_default, lc)
                          for D in D_grid])
    ax.plot(D_grid, W_C_curve, color=col, linewidth=2.0,
            label=f"$\\lambda_c = {lc:g}$", zorder=3)


ax.set_xlabel("Effective motility  $D$")
ax.set_ylabel("Contact-weapon fitness  $W_C(D)$")
ax.set_xlim(0, D_MAX_PLOT)
ax.set_title("Contact weapons (T6SS): robust gain", color=COL_T6SS,
             fontsize=10.0, pad=6)
ax.legend(loc="upper left", title="competitor\ndensity",
          title_fontsize=8.5)

panel_label(ax, "B")


# ---------------------------------------------------------------------------
# Panel C — Diffusible weapon W_D vs D, several lambda_s
# ---------------------------------------------------------------------------
ax = ax_C
ell_default = 0.40
ls_values = [0.1, 0.3, 0.6, 1.0]
# Use a perceptually distinct sequence: light to dark teal/blue
cmap = plt.get_cmap("viridis_r")
colors_C = [cmap(0.10), cmap(0.35), cmap(0.55), cmap(0.78)]

ax.axhline(0, color="#BBBBBB", linewidth=0.7, zorder=1)
for ls, col in zip(ls_values, colors_C):
    W_D_curve = np.array([diffusible_fitness(D, ls, ell_default)
                          for D in D_grid])
    ax.plot(D_grid, W_D_curve, color=col, linewidth=1.9,
            label=f"$\\lambda_s = {ls:g}$")


ax.set_xlabel("Effective motility  $D$")
ax.set_ylabel("Diffusible-weapon fitness  $W_D(D)$")
ax.set_xlim(0, D_MAX_PLOT)
ax.set_title("Diffusible weapons: context-dependent", color=COL_BACT,
             fontsize=10.0, pad=6)
ax.legend(loc="upper right", title=f"$\\ell = {ell_default:g}$",
          title_fontsize=8.5)

panel_label(ax, "C")


# ---------------------------------------------------------------------------
# Panel D — Direct vs kin decomposition of W_D
# ---------------------------------------------------------------------------
ax = ax_D
ls_for_decomp = 0.8
ell_for_decomp = 0.5
D_grid_D = np.linspace(0.001, D_MAX_PLOT, 280)
direct = np.zeros_like(D_grid_D)
kin    = np.zeros_like(D_grid_D)
total  = np.zeros_like(D_grid_D)
for i, D in enumerate(D_grid_D):
    d, k, t = diffusible_components(D, ls_for_decomp, ell_for_decomp)
    direct[i] = d
    kin[i] = k
    total[i] = t

ax.axhline(0, color="#BBBBBB", linewidth=0.7, zorder=1)
ax.fill_between(D_grid_D, 0, direct, color=COL_DIRECT, alpha=0.12, zorder=2)
ax.fill_between(D_grid_D, 0, kin,    color=COL_KIN,    alpha=0.12, zorder=2)
ax.plot(D_grid_D, direct, color=COL_DIRECT, linewidth=2.0,
        label="direct producer share", zorder=3)
ax.plot(D_grid_D, kin,    color=COL_KIN,    linewidth=2.0,
        label="kin-mediated share",    zorder=3)
ax.plot(D_grid_D, total,  color=COL_BACT,   linewidth=2.3, linestyle="--",
        label="net $W_D$ (minus cost)", zorder=4)

# Mark peak and place a label just below it in the empty space under the
# zero line. The label is centered horizontally on the peak.
ipk = int(np.argmax(total))
ax.scatter([D_grid_D[ipk]], [total[ipk]], color=COL_BACT, s=40,
           zorder=5, edgecolor="white", linewidth=0.7)
ax.set_ylim(0, max(kin.max(), total.max()) * 1.10)
# Place the peak label in the gap between the green (direct) curve and the
# dashed net curve, with a thin connector to the peak marker.
label_xy = (0.13, 0.093)
ax.annotate(
    f"net peak  $D \\approx {D_grid_D[ipk]:.2f}$",
    xy=(D_grid_D[ipk], total[ipk]),
    xytext=label_xy,
    fontsize=8.2, color=COL_BACT, ha="center", va="center",
    arrowprops=dict(arrowstyle="-", color=COL_BACT, linewidth=0.7,
                    shrinkA=2, shrinkB=4),
)


ax.set_xlabel("Effective motility  $D$")
ax.set_ylabel("Fitness contribution")
ax.set_xlim(0, D_MAX_PLOT)
ax.set_title("Why diffusible payoff peaks at intermediate motility",
             fontsize=10.0, pad=6)
ax.legend(loc="upper right",
          title=f"$\\lambda_s = {ls_for_decomp:g}$,  $\\ell = {ell_for_decomp:g}$",
          title_fontsize=8.5)

panel_label(ax, "D")


# ---------------------------------------------------------------------------
# Panel E — Heatmap of motility effect on diffusible weapons
# ---------------------------------------------------------------------------
ax = ax_E
ell_grid = np.linspace(0.05, 1.4, 90)
ls_grid  = np.linspace(0.02, 1.0, 90)
delta = np.zeros((len(ls_grid), len(ell_grid)))
for i, ls in enumerate(ls_grid):
    for j, ell in enumerate(ell_grid):
        delta[i, j] = (diffusible_fitness(D_m, ls, ell)
                       - diffusible_fitness(D_nm, ls, ell))

vmax = np.max(np.abs(delta))
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
im = ax.pcolormesh(ell_grid, ls_grid, delta, cmap="RdBu_r", norm=norm,
                   shading="auto")
# Zero contour
ax.contour(ell_grid, ls_grid, delta, levels=[0],
           colors="black", linewidths=0.9)

ax.set_xlabel("Toxin range  $\\ell$")
ax.set_ylabel("Conspecific density  $\\lambda_s$")
ax.set_title("Motility effect on diffusible weapons", fontsize=10.0, pad=6)

cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
cbar.set_label(r"$W_D(D_{\rm mot}) - W_D(D_{\rm non})$", fontsize=8.8)
cbar.ax.tick_params(labelsize=8)

# Sign labels — placed in the regions where the sign holds (red = motility
# helps, blue = motility hurts). White outline keeps them legible against
# saturated heatmap colours.
white_outline = [pe.withStroke(linewidth=2.4, foreground="white")]
ax.text(0.85, 0.10, "motility helps", color="#7A1F1F",
        fontsize=9.5, va="center", ha="center", fontweight="bold",
        path_effects=white_outline)
ax.text(1.30, 0.92, "motility hurts", color="#1F3F76",
        fontsize=9.5, va="top", ha="right", fontweight="bold",
        path_effects=white_outline)

panel_label(ax, "E")


# ---------------------------------------------------------------------------
# Panel F — Macroevolutionary contrast: distribution of motility-induced
# fitness gain across parameter space, T6SS vs bacteriocins
# ---------------------------------------------------------------------------
ax = ax_F

rng = np.random.default_rng(0)
n = 1500
# Independent log-normal sampling on each parameter (a multivariate Gaussian
# on log-parameters with diagonal covariance). Medians are set to biologically
# plausible central values and sigma=0.55 gives roughly a 3-fold 5-95 spread,
# placing most weight on typical lineages while still covering the range that
# panels C-E display.
median_ls, median_lc, median_ell = 0.30, 1.00, 0.40
sigma_log = 0.55
lambda_s_samples = np.exp(rng.normal(np.log(median_ls),  sigma_log, size=n))
lambda_c_samples = np.exp(rng.normal(np.log(median_lc),  sigma_log, size=n))
ell_samples      = np.exp(rng.normal(np.log(median_ell), sigma_log, size=n))
# Clip to the displayed parameter ranges so out-of-bound draws don't dominate
lambda_s_samples = np.clip(lambda_s_samples, 0.02, 1.5)
lambda_c_samples = np.clip(lambda_c_samples, 0.10, 4.0)
ell_samples      = np.clip(ell_samples,      0.05, 1.5)

deltaW_C = np.array([
    contact_fitness(D_m, ls, lc) - contact_fitness(D_nm, ls, lc)
    for ls, lc in zip(lambda_s_samples, lambda_c_samples)
])
deltaW_D = np.array([
    diffusible_fitness(D_m, ls, ell, lc) - diffusible_fitness(D_nm, ls, ell, lc)
    for ls, ell, lc in zip(lambda_s_samples, ell_samples, lambda_c_samples)
])

ax.axvline(0, color="#888888", linewidth=0.9, linestyle="--", zorder=1)

# Two horizontal strip+violin rows
groups = [
    ("contact (T6SS)",      deltaW_C, COL_T6SS, 1.0),
    ("diffusible\n(bacteriocin)", deltaW_D, COL_BACT, 0.0),
]

for name, vals, col, y in groups:
    # Violin
    parts = ax.violinplot([vals], positions=[y], vert=False,
                          widths=0.7, showmeans=False, showmedians=False,
                          showextrema=False)
    for body in parts["bodies"]:
        body.set_facecolor(col)
        body.set_alpha(0.30)
        body.set_edgecolor(col)
        body.set_linewidth(0.8)

    # Jittered strip
    y_jit = y + rng.uniform(-0.13, 0.13, size=len(vals))
    ax.scatter(vals, y_jit, s=6, color=col, alpha=0.30,
               edgecolor="none", zorder=2)

    # Median + IQR
    q25, q50, q75 = np.percentile(vals, [25, 50, 75])
    ax.plot([q25, q75], [y, y], color=col, linewidth=2.6, zorder=4)
    ax.scatter([q50], [y], color="white", edgecolor=col,
               s=42, linewidth=1.5, zorder=5)

# Annotate median and P(>0) above each violin. Both labels are aligned on
# the same x position (the contact T6SS median) so they read as a column.
contact_med = np.median(deltaW_C)
for name, vals, col, y in groups:
    frac = (vals > 0).mean() * 100
    med = np.median(vals)
    label = f"median = {med:+.2f},   $P(\\Delta W>0) = {frac:.0f}\\%$"
    ax.text(contact_med, y + 0.46, label,
            fontsize=8.4, color=col, ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      edgecolor="none"))

ax.set_yticks([0, 1])
ax.set_yticklabels(["diffusible\n(bacteriocin)", "contact\n(T6SS)"])
ax.set_ylim(-0.5, 1.7)
# Widen x range so the bacteriocin label (centered near 0) fits inside.
xpad = 0.30
xmin = min(deltaW_C.min(), deltaW_D.min()) - xpad
xmax = max(deltaW_C.max(), deltaW_D.max()) + 0.10
ax.set_xlim(xmin, xmax)
ax.set_xlabel(r"Motility-induced fitness gain  $\Delta W = W(D_{\rm mot}) - W(D_{\rm non})$")
ax.set_title("Predicted macroevolutionary signal", fontsize=10.0, pad=6)

panel_label(ax, "F")


# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out_pdf = os.path.join(HERE, "model_figure.pdf")
out_png = os.path.join(HERE, "model_figure.png")
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
fig.savefig(out_png, dpi=300, bbox_inches="tight")
print(f"wrote {out_pdf}")
print(f"wrote {out_png}")
