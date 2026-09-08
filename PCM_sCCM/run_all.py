"""
run_all.py -- driver script: runs Model A (0-D lumped sCCM force-balance)
and Model B (1-D radial transient enthalpy PDE) across the three named
slip-length scenarios in parameters.py, prints a milestone summary table,
and saves comparison plots to results/.

Usage:
    python3 run_all.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from parameters import Params
import ccm_slip_model as modelA
import radial_enthalpy_model as modelB

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

COLORS = {"no slip (b=0, baseline CCM)": "#4C6EF5",
          "slip b=45 um (paper, low end)": "#F59F00",
          "slip b=90 um (paper, high end)": "#E03131"}


def fmt_t(x):
    if x is None:
        return "not reached"
    return f"{x:6.1f} s ({x/60.0:5.2f} min)"


def main():
    P = Params()

    print("=" * 78)
    print("sCCM vertical-tube pulsed-heater model -- Model A (0-D lumped CCM)")
    print("=" * 78)
    resultsA = {}
    for name, b in P.slip_scenarios.items():
        resultsA[name] = modelA.simulate(P, b, verbose=False)

    print(f"{'scenario':38s} {'t50':>18s} {'t90':>18s} {'t95':>18s} {'t99':>18s}")
    for name, res in resultsA.items():
        m = res.milestones
        print(f"{name:38s} {fmt_t(m['t50']):>18s} {fmt_t(m['t90']):>18s} {fmt_t(m['t95']):>18s} {fmt_t(m['t99']):>18s}")

    baseline_t90 = resultsA["no slip (b=0, baseline CCM)"].milestones["t90"]
    print()
    for name, res in resultsA.items():
        t90 = res.milestones["t90"]
        if baseline_t90 and t90:
            print(f"  {name}: t90 speed-up vs no-slip = {baseline_t90/t90:5.2f}x")

    print()
    print("=" * 78)
    print("sCCM vertical-tube pulsed-heater model -- Model B (1-D radial enthalpy PDE)")
    print("=" * 78)
    resultsB = {}
    for name, b in P.slip_scenarios.items():
        resultsB[name] = modelB.simulate(P, b, verbose=False)

    print(f"{'scenario':38s} {'t50':>18s} {'t90':>18s} {'t95':>18s} {'t99':>18s}")
    for name, res in resultsB.items():
        m = res.milestones
        print(f"{name:38s} {fmt_t(m['t50']):>18s} {fmt_t(m['t90']):>18s} {fmt_t(m['t95']):>18s} {fmt_t(m['t99']):>18s}")

    # ---- Plot 1: Model A melt fraction vs time -------------------------
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for name, res in resultsA.items():
        ax.plot(res.t / 60.0, 100 * res.melt_fraction, label=name, color=COLORS[name], lw=1.8)
    ax.set_xlabel("time [min]")
    ax.set_ylabel("PCM melted [%] (Model A, sinking-column mass basis)")
    ax.set_title("Model A: sCCM lumped force-balance -- melt progress")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, "A_melt_fraction_vs_time.png"), dpi=150)
    plt.close(fig)

    # ---- Plot 2: Model A film thickness & heat flux (first 10 min) ----
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 7), sharex=True)
    tcut = 600.0
    for name, res in resultsA.items():
        mask = res.t <= tcut
        ax1.plot(res.t[mask], res.delta[mask] * 1e6, label=name, color=COLORS[name], lw=1.4)
        ax2.plot(res.t[mask], res.q_h[mask] / 1000.0, label=name, color=COLORS[name], lw=1.4)
    ax1.set_ylabel("film thickness delta [um]")
    ax1.set_title("Model A: CCM film thickness and heater-surface heat flux (first 10 min)")
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=8)
    ax2.set_ylabel("heat flux q_h [kW/m^2]")
    ax2.set_xlabel("time [s]")
    ax2.grid(alpha=0.3)
    # shade OFF intervals using the baseline scenario's pulse trace
    base = resultsA["no slip (b=0, baseline CCM)"]
    mask = base.t <= tcut
    off = ~base.pulse_on[mask]
    ax2.fill_between(base.t[mask], 0, ax2.get_ylim()[1], where=off, color="grey", alpha=0.15, step="mid")
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, "A_film_and_flux.png"), dpi=150)
    plt.close(fig)

    # ---- Plot 3: Model B liquid fraction vs time -----------------------
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for name, res in resultsB.items():
        ax.plot(res.t / 60.0, 100 * res.liquid_frac_vol, label=name, color=COLORS[name], lw=1.8)
    ax.set_xlabel("time [min]")
    ax.set_ylabel("volume-averaged liquid fraction [%] (Model B, fixed annulus)")
    ax.set_title("Model B: 1-D radial enthalpy PDE -- melt progress")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, "B_liquid_fraction_vs_time.png"), dpi=150)
    plt.close(fig)

    # ---- Plot 4: Model B melt front vs time, pulsing visible -----------
    fig, ax = plt.subplots(figsize=(7.5, 5))
    tcut_B = 300.0
    for name, res in resultsB.items():
        mask = res.t <= tcut_B
        ax.plot(res.t[mask], (res.r_melt[mask] - P.r_h) * 1000.0, label=name, color=COLORS[name], lw=1.6)
    base = resultsB["no slip (b=0, baseline CCM)"]
    mask = base.t <= tcut_B
    off = ~base.pulse_on[mask]
    ax.fill_between(base.t[mask], 0, ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 1,
                     where=off, color="grey", alpha=0.15, step="mid", label="pulse OFF")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("melt-layer thickness  r_melt - r_h  [mm]")
    ax.set_title("Model B: melt-front growth, pulsing stalls visible (first 5 min)")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, "B_melt_front_pulsing.png"), dpi=150)
    plt.close(fig)

    # ---- Plot 5: Model B temperature-field snapshots (default slip) ---
    default_name = list(P.slip_scenarios.keys())[1]  # the paper's low-end slip scenario, b=45 um
    res = resultsB[default_name]
    fig, ax = plt.subplots(figsize=(7.5, 5))
    n_snap = res.T_history.shape[0]
    idxs = np.linspace(0, n_snap - 1, min(8, n_snap)).astype(int)
    cmap = plt.get_cmap("viridis")
    for j, idx in enumerate(idxs):
        ax.plot((res.r_c - P.r_h) * 1000.0, res.T_history[idx], color=cmap(j / max(1, len(idxs) - 1)),
                label=f"t={res.t_saved[idx]:.0f} s")
    ax.axhline(P.Ts, color="k", ls=":", lw=0.8, label="Ts, Tl (mushy band)")
    ax.axhline(P.Tl, color="k", ls=":", lw=0.8)
    ax.set_xlabel("radial distance from heater surface  r - r_h  [mm]")
    ax.set_ylabel("temperature [deg C]")
    ax.set_title(f"Model B: radial temperature profiles ({default_name})")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, "B_temperature_profiles.png"), dpi=150)
    plt.close(fig)

    print()
    print(f"Plots written to {RESULTS_DIR}/")


if __name__ == "__main__":
    main()
