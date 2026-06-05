"""
The sensitivity-analysis figures for Supplementary Methods S1.

  Fig S1 — Spatial scales L(D) and xi(D) as functions of effective motility.
           Shows the two competing spatial effects of motility as continuous
           curves and as a relatedness-vs-distance plot.
  Fig S2 — Inclusive-fitness value G(r, D, lambda_s) of a kill as a function
           of kill distance, for low vs high motility and low vs high
           conspecific density.
  Fig S3 — Robustness of the sign-flip surface in (ell, lambda_s) across a
           range of competitor densities lambda_c.
  Fig S4 — Sensitivity of the diffusible-weapon motility effect to the
           relatedness parameters xi_0 (clonal patch length) and eta
           (motility-driven kin erosion).
"""

import os
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from kin_competition_motility_weapon_model import (
    params,
    exploitation_length,
    relatedness_length,
    relatedness_at_distance,
    conspecific_competition_weight,
    focal_access_weight,
    per_kill_inclusive_benefit,
    diffusible_kill_distance_pdf,
    diffusible_total_kill_opportunity,
    diffusible_fitness,
)


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9.5,
    "axes.labelsize": 10,
    "axes.titlesize": 10.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.9,
    "legend.frameon": False,
    "legend.fontsize": 8.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

COL_NONMOT = "#2E5E8C"
COL_MOT    = "#D7642C"
COL_L      = "#1F7A4C"   # exploitation length L
COL_XI     = "#7B3294"   # relatedness scale xi

D_nm = params["D_nonmotile"]
D_m  = params["D_motile"]
D_MAX = 0.7


def panel_label(ax, label, x=-0.18, y=1.06):
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="bottom", ha="left")


# ---------------------------------------------------------------------------
# Fig S2 — Spatial scales L(D) and xi(D)
# ---------------------------------------------------------------------------
def fig_s2():
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6))

    # Panel A: L(D) and xi(D) curves
    ax = axes[0]
    D_grid = np.linspace(0.001, D_MAX, 240)
    L = np.array([exploitation_length(D) for D in D_grid])
    xi = np.array([relatedness_length(D) for D in D_grid])
    ax.plot(D_grid, L,  color=COL_L,  linewidth=2.2, label="exploitation reach $L(D)$")
    ax.plot(D_grid, xi, color=COL_XI, linewidth=2.2, label="kin scale $\\xi(D)$")
    ax.scatter([D_nm, D_m], [exploitation_length(D_nm), exploitation_length(D_m)],
               color=COL_L, s=36, zorder=4, edgecolor="white", linewidth=0.7)
    ax.scatter([D_nm, D_m], [relatedness_length(D_nm), relatedness_length(D_m)],
               color=COL_XI, s=36, zorder=4, edgecolor="white", linewidth=0.7)
    ax.set_xlabel("Effective motility  $D$")
    ax.set_ylabel("Spatial scale")
    ax.set_xlim(0, D_MAX)
    ax.set_ylim(0, max(L.max(), xi.max()) * 1.1)
    ax.set_title("Two spatial scales set by motility", fontsize=10, pad=6)
    ax.legend(loc="upper right")
    panel_label(ax, "A")

    # Panel B: R(r, D) for non-motile and motile
    ax = axes[1]
    r_grid = np.linspace(0, 2.5, 300)
    R_nm = relatedness_at_distance(r_grid, D_nm)
    R_m  = relatedness_at_distance(r_grid, D_m)
    ax.plot(r_grid, R_nm, color=COL_NONMOT, linewidth=2.2,
            label=f"non-motile ($D = {D_nm:g}$)")
    ax.plot(r_grid, R_m,  color=COL_MOT,    linewidth=2.2,
            label=f"motile ($D = {D_m:g}$)")
    ax.axhline(0, color="#BBB", linewidth=0.7)
    ax.set_xlabel("Distance from focal producer  $r$")
    ax.set_ylabel("Relatedness  $R(r, D)$")
    ax.set_xlim(0, 2.5)
    ax.set_ylim(0, params["max_relatedness"] * 1.05)
    ax.set_title("Relatedness decays faster under high motility",
                 fontsize=10, pad=6)
    ax.legend(loc="upper right")
    panel_label(ax, "B")

    plt.tight_layout()
    out = os.path.join(HERE, "fig_S2_spatial_scales.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(out.replace(".png", ".pdf"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Fig S3 — Inclusive value G(r, D, lambda_s) vs kill distance
# ---------------------------------------------------------------------------
def fig_s3():
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6))

    r_grid = np.linspace(0.001, 2.5, 400)

    # Panel A: low conspecific density
    ax = axes[0]
    ls = 0.3
    G_nm = np.array([per_kill_inclusive_benefit(r, D_nm, ls) for r in r_grid])
    G_m  = np.array([per_kill_inclusive_benefit(r, D_m,  ls) for r in r_grid])
    ax.plot(r_grid, G_nm, color=COL_NONMOT, linewidth=2.2,
            label=f"non-motile ($D = {D_nm:g}$)")
    ax.plot(r_grid, G_m,  color=COL_MOT,    linewidth=2.2,
            label=f"motile ($D = {D_m:g}$)")
    ax.set_xlabel("Kill distance  $r$")
    ax.set_ylabel("Inclusive value  $G(r, D, \\lambda_s)$")
    ax.set_xlim(0, 2.5)
    ax.set_ylim(0, 1.05)
    ax.set_title(f"Sparse conspecifics  ($\\lambda_s = {ls:g}$)",
                 fontsize=10, pad=6)
    ax.legend(loc="upper right")
    panel_label(ax, "A")

    # Panel B: high conspecific density
    ax = axes[1]
    ls = 1.0
    G_nm = np.array([per_kill_inclusive_benefit(r, D_nm, ls) for r in r_grid])
    G_m  = np.array([per_kill_inclusive_benefit(r, D_m,  ls) for r in r_grid])
    ax.plot(r_grid, G_nm, color=COL_NONMOT, linewidth=2.2,
            label=f"non-motile ($D = {D_nm:g}$)")
    ax.plot(r_grid, G_m,  color=COL_MOT,    linewidth=2.2,
            label=f"motile ($D = {D_m:g}$)")
    ax.set_xlabel("Kill distance  $r$")
    ax.set_ylabel("Inclusive value  $G(r, D, \\lambda_s)$")
    ax.set_xlim(0, 2.5)
    ax.set_ylim(0, 1.05)
    ax.set_title(f"Dense conspecifics  ($\\lambda_s = {ls:g}$)",
                 fontsize=10, pad=6)
    ax.legend(loc="upper right")
    panel_label(ax, "B")

    plt.tight_layout()
    out = os.path.join(HERE, "fig_S3_inclusive_value.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(out.replace(".png", ".pdf"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Fig S4 — Robust sign-flip across competitor density lambda_c
# ---------------------------------------------------------------------------
def fig_s4():
    lc_values = [0.3, 1.0, 3.0]
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.8))

    ell_grid = np.linspace(0.05, 1.4, 70)
    ls_grid  = np.linspace(0.02, 1.0, 70)

    # First compute all to get common color scale
    grids = []
    for lc in lc_values:
        delta = np.zeros((len(ls_grid), len(ell_grid)))
        for i, ls in enumerate(ls_grid):
            for j, ell in enumerate(ell_grid):
                delta[i, j] = (
                    diffusible_fitness(D_m, ls, ell, lc)
                    - diffusible_fitness(D_nm, ls, ell, lc)
                )
        grids.append(delta)
    vmax = max(np.max(np.abs(g)) for g in grids)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    for ax, lc, delta, lbl in zip(axes, lc_values, grids, ["A", "B", "C"]):
        im = ax.pcolormesh(ell_grid, ls_grid, delta, cmap="RdBu_r",
                           norm=norm, shading="auto")
        ax.contour(ell_grid, ls_grid, delta, levels=[0],
                   colors="black", linewidths=0.8)
        ax.set_xlabel("Toxin range  $\\ell$")
        ax.set_ylabel("Conspecific density  $\\lambda_s$")
        ax.set_title(f"$\\lambda_c = {lc:g}$", fontsize=10, pad=6)
        panel_label(ax, lbl)

    cbar = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label(r"$W_D(D_{\rm mot}) - W_D(D_{\rm non})$", fontsize=8.8)
    cbar.ax.tick_params(labelsize=8)

    out = os.path.join(HERE, "fig_S4_lambda_c_robust.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(out.replace(".png", ".pdf"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Fig S5 — Sensitivity to relatedness parameters xi_0 and eta
# ---------------------------------------------------------------------------
def fig_s5():
    # Vary xi_0 and eta around the baseline used in Fig. 1.
    baseline_xi0 = params["nonmotile_relatedness_length"]   # 1.20
    baseline_eta = params["motility_relatedness_decay"]     # 4.00

    # Sweep ell and lambda_s, compute delta W_D for several relatedness configs.
    ell_grid = np.linspace(0.05, 1.4, 60)
    ls_grid  = np.linspace(0.02, 1.0, 60)
    lc = 1.0

    configs = [
        ("low kin structure",  baseline_xi0 * 0.5, baseline_eta * 2.0),
        ("baseline",           baseline_xi0,       baseline_eta),
        ("high kin structure", baseline_xi0 * 1.5, baseline_eta * 0.5),
    ]

    def delta_with_relatedness(xi0, eta):
        # Temporarily override params, restore after.
        old_xi0 = params["nonmotile_relatedness_length"]
        old_eta = params["motility_relatedness_decay"]
        params["nonmotile_relatedness_length"] = xi0
        params["motility_relatedness_decay"]   = eta
        delta = np.zeros((len(ls_grid), len(ell_grid)))
        for i, ls in enumerate(ls_grid):
            for j, ell in enumerate(ell_grid):
                delta[i, j] = (
                    diffusible_fitness(D_m, ls, ell, lc)
                    - diffusible_fitness(D_nm, ls, ell, lc)
                )
        params["nonmotile_relatedness_length"] = old_xi0
        params["motility_relatedness_decay"]   = old_eta
        return delta

    grids = [delta_with_relatedness(xi0, eta) for _, xi0, eta in configs]
    vmax = max(np.max(np.abs(g)) for g in grids)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.8))
    for ax, (title, xi0, eta), delta, lbl in zip(
        axes, configs, grids, ["A", "B", "C"]
    ):
        im = ax.pcolormesh(ell_grid, ls_grid, delta, cmap="RdBu_r",
                           norm=norm, shading="auto")
        ax.contour(ell_grid, ls_grid, delta, levels=[0],
                   colors="black", linewidths=0.8)
        ax.set_xlabel("Toxin range  $\\ell$")
        ax.set_ylabel("Conspecific density  $\\lambda_s$")
        ax.set_title(f"{title}\n$\\xi_0 = {xi0:.2f}$, $\\eta = {eta:.2f}$",
                     fontsize=9.5, pad=6)
        panel_label(ax, lbl)

    cbar = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label(r"$W_D(D_{\rm mot}) - W_D(D_{\rm non})$", fontsize=8.8)
    cbar.ax.tick_params(labelsize=8)

    out = os.path.join(HERE, "fig_S5_relatedness_sensitivity.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(out.replace(".png", ".pdf"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


if __name__ == "__main__":
    for fn in [fig_s2, fig_s3, fig_s4, fig_s5]:
        out = fn()
        print(f"wrote {out}")
