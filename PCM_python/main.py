"""
main.py -- top-level driver for the "hot pipe / PCM / pulsed copper
sleeve / PCM / insulated wall" case (see PCM_python/README.md).

Runs:
  1) the case exactly as drawn (copper sleeve centered in the pipe->wall
     gap), with the sleeve pulse-heated at the lowest peak power that
     reaches the PCM melt point within the target ramp time;
  2) a no-sleeve (pipe-only) reference, for comparison;
  3) a sweep of the sleeve's radial position (same "lowest power"
     control rule at every position) to find the minimum achievable
     melting time for this architecture;
and writes plots + a CSV/markdown summary to PCM_python/results/.
"""

import csv
import json
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from params import Params
from model import Grid
from simulate import build_case, find_min_power, run_full_melt
from sweep import run_position_sweep, run_no_sleeve_reference

RESULTS = "results"


def hms(sec):
    if sec is None or (isinstance(sec, float) and np.isnan(sec)):
        return "n/a"
    h = int(sec // 3600)
    m = int((sec % 3600) // 60)
    s = int(sec % 60)
    return f"{h:d}h {m:02d}m {s:02d}s ({sec:.0f} s)"


def main():
    P = Params()
    print("=== PCM radial model: hot pipe + pulsed copper sleeve ===")
    print(f"r_pipe={P.r_pipe*1000:.2f} mm, r_wall={P.r_wall*1000:.2f} mm "
          f"(gap width {(P.r_wall-P.r_pipe)*1000:.2f} mm) [container size ASSUMED]")
    print(f"T_hotwater={P.T_hotwater} C (Dirichlet at pipe surface), Tm={P.Tm} C, "
          f"T_init={P.T_init} C")
    print(f"PCM: rho={P.rho_pcm} kg/m3, cp={P.cp_pcm} J/kgK, k_s={P.k_pcm_s}, "
          f"k_l={P.k_pcm_l} W/mK, L={P.L_pcm} J/kg  [material dataset ASSUMED -- "
          f"paraffin/RT42-class placeholder]")
    print(f"Copper sleeve thickness={P.t_sleeve*1000:.2f} mm, ramp target="
          f"{P.t_ramp_target:.0f} s\n")

    # ---- 1) Baseline: sleeve centered in the gap, as drawn -----------------
    print("--- Baseline: sleeve centered in the gap (as drawn) ---")
    grid, is_cu, (r_in, r_out) = build_case(P, center_frac=0.5)
    t0 = time.time()
    P_min_base, T_reached = find_min_power(P, grid, is_cu)
    res_base = run_full_melt(P, grid, is_cu, P_min_base, record_every=20, store_snapshots=True)
    print(f"  sleeve band: {r_in*1000:.2f}-{r_out*1000:.2f} mm")
    print(f"  minimum sleeve peak power P_min = {P_min_base:.1f} W/m of pipe "
          f"(reaches {T_reached:.1f} C by t_ramp={P.t_ramp_target:.0f} s)")
    m = res_base["milestones"]
    print(f"  t50={hms(m['t50'])}  t90={hms(m['t90'])}  t99={hms(m['t99'])}  "
          f"t999(full melt)={hms(m['t999'])}")
    print(f"  sleeve heater duty cycle over the whole run: {res_base['duty_cycle']*100:.1f}%")
    print(f"  elapsed wall time: {time.time()-t0:.1f} s\n")

    # ---- 2) No-sleeve reference --------------------------------------------
    print("--- Reference: no sleeve, pipe conduction only ---")
    _, res_ref = run_no_sleeve_reference(P, record_every=20)
    mr = res_ref["milestones"]
    print(f"  t50={hms(mr['t50'])}  t90={hms(mr['t90'])}  t99={hms(mr['t99'])}  "
          f"t999(full melt)={hms(mr['t999'])}\n")

    # ---- 3) Sleeve-position sweep -------------------------------------------
    print("--- Sleeve-position sweep (coarse) ---")
    fracs_coarse = np.linspace(0.15, 0.90, 16)
    rows = run_position_sweep(P, fracs=fracs_coarse, record_every=40)

    best_coarse = min(rows, key=lambda r: r["t999"] if r["t999"] is not None else 1e18)
    print(f"\n  coarse best: frac={best_coarse['frac']:.2f} -> t999={best_coarse['t999']:.0f} s")

    print("\n--- Sleeve-position sweep (refine near coarse optimum) ---")
    lo = max(0.10, best_coarse["frac"] - 0.08)
    hi = min(0.95, best_coarse["frac"] + 0.08)
    fracs_fine = np.linspace(lo, hi, 9)
    rows_fine = run_position_sweep(P, fracs=fracs_fine, record_every=40)
    rows_all = rows + rows_fine
    best = min(rows_all, key=lambda r: r["t999"] if r["t999"] is not None else 1e18)

    print(f"\n  OPTIMAL sleeve position: frac={best['frac']:.3f} of the pipe->wall gap "
          f"(r_mid={best['r_mid']*1000:.2f} mm)")
    print(f"  P_min at optimum = {best['P_min_W_per_m']:.1f} W/m")
    print(f"  t99={hms(best['t99'])}   t999(full melt)={hms(best['t999'])}")

    # ---- Full run at the optimal position (for plotting) --------------------
    grid_opt, is_cu_opt, (r_in_opt, r_out_opt) = build_case(P, center_frac=best["frac"])
    res_opt = run_full_melt(P, grid_opt, is_cu_opt, best["P_min_W_per_m"],
                             record_every=20, store_snapshots=True)

    # ---- Save CSV of the sweep ------------------------------------------------
    import os
    os.makedirs(RESULTS, exist_ok=True)
    with open(f"{RESULTS}/sweep.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_all[0].keys()))
        w.writeheader()
        for r in sorted(rows_all, key=lambda r: r["frac"]):
            w.writerow(r)

    # ---- Plots ------------------------------------------------------------
    make_plots(P, res_base, res_ref, res_opt, rows_all, best,
               (r_in, r_out), (r_in_opt, r_out_opt))

    # ---- Summary ------------------------------------------------------------
    summary = dict(
        params=dict(r_pipe_mm=P.r_pipe*1000, r_wall_mm=P.r_wall*1000,
                     t_sleeve_mm=P.t_sleeve*1000, T_hotwater=P.T_hotwater,
                     Tm=P.Tm, T_init=P.T_init, t_ramp_target=P.t_ramp_target),
        baseline_centered=dict(P_min_W_per_m=P_min_base, milestones=m,
                                duty_cycle=res_base["duty_cycle"]),
        no_sleeve_reference=dict(milestones=mr),
        optimal=dict(frac=best["frac"], r_mid_mm=best["r_mid"]*1000,
                     P_min_W_per_m=best["P_min_W_per_m"],
                     t99=best["t99"], t999=best["t999"]),
    )
    with open(f"{RESULTS}/summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)

    print(f"\nWrote plots, sweep.csv and summary.json to {RESULTS}/")
    print("\n=== HEADLINE RESULT ===")
    print(f"As drawn (sleeve centered, {P_min_base:.0f} W/m minimum pulse power): "
          f"full melt in {hms(m['t999'])}")
    print(f"No auxiliary sleeve heater at all (pipe conduction only): "
          f"full melt in {hms(mr['t999'])}")
    print(f"Minimum achievable melting time for this architecture "
          f"(sleeve at {best['frac']*100:.0f}% of the way from pipe to wall, "
          f"{best['P_min_W_per_m']:.0f} W/m): full melt in {hms(best['t999'])}")


def make_plots(P, res_base, res_ref, res_opt, rows_all, best, band_base, band_opt):
    import os
    os.makedirs(RESULTS, exist_ok=True)

    # 1) Liquid fraction history: baseline vs optimal vs no-sleeve
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for res, label, style in (
        (res_base, "Sleeve centered in gap (as drawn)", "-"),
        (res_opt, f"Optimal position ({best['frac']*100:.0f}% pipe→wall)", "-"),
        (res_ref, "No sleeve (pipe conduction only)", "--"),
    ):
        h = res["hist"]
        ax.plot(h["t"] / 60.0, h["favg"] * 100, style, label=label, linewidth=1.8)
    ax.set_xlabel("time [min]")
    ax.set_ylabel("PCM average liquid fraction [%]")
    ax.set_title("Melting progress: effect of the pulsed copper sleeve")
    ax.axhline(99, color="gray", linewidth=0.7, linestyle=":")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{RESULTS}/fig1_liquid_fraction.png", dpi=150)
    plt.close(fig)

    # 2) Sweep curve: t99/t999 vs sleeve position
    rows_sorted = sorted(rows_all, key=lambda r: r["frac"])
    fracs = [r["frac"] for r in rows_sorted]
    t99s = [r["t99"] / 60.0 if r["t99"] else np.nan for r in rows_sorted]
    t999s = [r["t999"] / 60.0 if r["t999"] else np.nan for r in rows_sorted]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(fracs, t99s, "o-", label="t99 (99% melted)", markersize=4)
    ax.plot(fracs, t999s, "s-", label="t99.9 (fully melted)", markersize=4)
    ax.axvline(best["frac"], color="green", linestyle=":", linewidth=1.2,
               label=f"optimum ({best['frac']:.2f})")
    ax.axvline(0.5, color="gray", linestyle=":", linewidth=1.2, label="centered (as drawn)")
    ax.set_xlabel("sleeve position: fraction of pipe→wall gap")
    ax.set_ylabel("melting time [min]")
    ax.set_title("Melting time vs. copper-sleeve radial position")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{RESULTS}/fig2_position_sweep.png", dpi=150)
    plt.close(fig)

    # 3) Temperature profile snapshots (optimal case)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    grid_opt, is_cu_opt, _ = build_case(P, center_frac=best["frac"])
    snaps = res_opt["snapshots"]
    n_snap = len(snaps)
    pick = np.linspace(0, n_snap - 1, min(7, n_snap)).astype(int)
    cmap = plt.cm.viridis(np.linspace(0, 1, len(pick)))
    for c, i in zip(cmap, pick):
        s = snaps[i]
        ax.plot(grid_opt.r_c * 1000, s["T"], color=c, linewidth=1.4,
                label=f"t={s['t']/60:.0f} min")
    ax.axhline(P.Tm, color="red", linestyle="--", linewidth=1, label="PCM melt point")
    r_in, r_out = band_opt
    ax.axvspan(r_in * 1000, r_out * 1000, color="orange", alpha=0.25, label="copper sleeve")
    ax.set_xlabel("radius [mm]")
    ax.set_ylabel("temperature [°C]")
    ax.set_title("Radial temperature profiles (optimal sleeve position)")
    ax.legend(fontsize=7, loc="lower right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{RESULTS}/fig3_temperature_profiles.png", dpi=150)
    plt.close(fig)

    # 4) Sleeve temperature + pulsed power (baseline, centered)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 5.5), sharex=True)
    h = res_base["hist"]
    ax1.plot(h["t"] / 60.0, h["T_sleeve"], color="darkorange", linewidth=1.4)
    ax1.axhline(P.Tm, color="red", linestyle="--", linewidth=1)
    ax1.set_ylabel("sleeve avg. temp [°C]")
    ax1.grid(alpha=0.3)
    ax1.set_title("Sleeve pulse-heating behaviour (centered case)")
    ax2.plot(h["t"] / 60.0, h["power"], color="steelblue", linewidth=1.2, drawstyle="steps-post")
    ax2.set_xlabel("time [min]")
    ax2.set_ylabel("sleeve power [W/m]")
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{RESULTS}/fig4_sleeve_pulsing.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
