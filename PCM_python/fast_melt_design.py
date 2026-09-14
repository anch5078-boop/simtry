"""
fast_melt_design.py -- optimize the annulus dimensions (holding the
pipe radius and PCM material fixed) so the whole PCM charge melts in
under 3 minutes, using the same pipe + PCM + pulsed copper sleeve +
insulated wall architecture as main.py, but re-targeted at melt SPEED
rather than minimum standby power.

Key difference from main.py's case: there, the sleeve's job was to sit
at the melt point using the lowest power that could get it there
("efficiency" framing). Here the objective changed to "melt everything
in <180 s", so the sleeve is run near its full available power for the
whole transient (T_setpoint below caps it only as a safety ceiling,
not an efficiency target -- see the module docstring debate in
simulate.run_full_melt).

Why dimensions have to shrink a lot: conduction time through a low-
conductivity PCM (k ~ 0.2 W/m*K) scales with the SQUARE of the distance
heat has to travel. The original 32.3 mm pipe->wall gap took ~2-3.5 h
to melt; getting under 180 s is not a tweak, it is roughly a (180 /
8500)^0.5 ~ 1/7 reduction in the relevant length scale -- consistent
with the sweep below, where even pure pipe conduction alone (no
sleeve) needs the gap down to ~4 mm to melt in <180 s.

This script:
  1) sweeps gap width (pipe->wall spacing) for the pure-conduction
     (no sleeve) case, to find the "no auxiliary heater needed" limit;
  2) sweeps gap width x sleeve peak power (sleeve position fixed near
     its own local optimum, found by a quick position sweep) to find,
     for each gap, the minimum sleeve power that still clears the
     <180 s target -- this quantifies how much extra PCM thickness
     (= extra thermal storage) the pulsed sleeve buys, for how much
     power;
  3) locks in one recommended design (documented, with margin under
     180 s and a safety cap on sleeve temperature) and produces the
     detailed plots for it.
"""

import csv
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from params import Params
from model import Grid
from simulate import build_case, run_full_melt

RESULTS = "results_fast_melt"
TARGET_S = 180.0          # the <3-minute requirement
BISECT_TARGET_S = 170.0   # bisect to a bit inside the target for margin
T_SAFETY_CAP = 180.0      # deg C -- sleeve heater cuts off above this
                          # (well below paraffin thermal-degradation
                          # territory, ~200-250 C) -- ASSUMED safety bound
GAP_STEPS_NO_SLEEVE = [2, 3, 4, 5, 6, 8, 10, 12]
GAP_STEPS_SLEEVE = [6, 7, 8, 9, 10, 11, 12, 13, 14]


def t999_of(P, grid, is_cu, power, record_every=100):
    res = run_full_melt(P, grid, is_cu, power, record_every=record_every, T_setpoint=T_SAFETY_CAP)
    return res["milestones"]["t999"], res


def sweep_no_sleeve():
    rows = []
    for gap_mm in GAP_STEPS_NO_SLEEVE:
        P = Params(r_wall=0.0127 + gap_mm / 1000.0, n_cells=220, dt=0.25, t_max=1200)
        grid = Grid(P.r_pipe, P.r_wall, P.n_cells)
        is_cu = np.zeros(grid.n, dtype=bool)
        t999, res = t999_of(P, grid, is_cu, 0.0, record_every=40)
        rows.append(dict(gap_mm=gap_mm, t999=t999))
        print(f"  [no sleeve] gap={gap_mm:3d} mm  t999={t999}")
    return rows


def position_sweep(P, fracs=None):
    if fracs is None:
        fracs = np.linspace(0.3, 0.85, 12)
    best_frac, best_t999 = None, np.inf
    for frac in fracs:
        grid, is_cu, _ = build_case(P, center_frac=frac)
        t999, _ = t999_of(P, grid, is_cu, 1500.0)  # fixed reference power just to rank positions
        if t999 is not None and t999 < best_t999:
            best_t999, best_frac = t999, frac
    return best_frac


def min_power_for_target(P, grid, is_cu, target=BISECT_TARGET_S, p_hi0=50.0, p_cap=20000.0, max_iter=28):
    p_hi = p_hi0
    while True:
        t999, _ = t999_of(P, grid, is_cu, p_hi)
        if t999 is not None and t999 <= target:
            break
        p_hi *= 1.6
        if p_hi > p_cap:
            return None, None
    p_lo = 0.0
    for _ in range(max_iter):
        mid = 0.5 * (p_lo + p_hi)
        t999, _ = t999_of(P, grid, is_cu, mid)
        if t999 is not None and t999 <= target:
            p_hi = mid
        else:
            p_lo = mid
    t999, _ = t999_of(P, grid, is_cu, p_hi)
    return p_hi, t999


def sweep_with_sleeve():
    rows = []
    for gap_mm in GAP_STEPS_SLEEVE:
        P = Params(r_wall=0.0127 + gap_mm / 1000.0, n_cells=260, dt=0.2, t_max=600, t_sleeve=0.001)
        frac = position_sweep(P)
        grid, is_cu, band = build_case(P, center_frac=frac)
        p_min, t999 = min_power_for_target(P, grid, is_cu)
        rows.append(dict(gap_mm=gap_mm, frac=frac, P_min_W_per_m=p_min, t999=t999))
        print(f"  [with sleeve] gap={gap_mm:3d} mm  frac={frac:.2f}  "
              f"P_min={p_min}  t999={t999}")
    return rows


def main():
    os.makedirs(RESULTS, exist_ok=True)

    print("=== Step 1: pure pipe conduction (no sleeve) -- gap-width sweep ===")
    rows_ns = sweep_no_sleeve()

    print("\n=== Step 2: with pulsed copper sleeve -- gap-width sweep "
          f"(bisecting sleeve power for t999 <= {BISECT_TARGET_S:.0f} s, "
          f"sleeve capped at {T_SAFETY_CAP:.0f} C) ===")
    rows_s = sweep_with_sleeve()

    # ---- Pick the recommended design: largest gap in the sleeve sweep whose
    # required power is still a modest, realistic value. -----------------
    feasible = [r for r in rows_s if r["P_min_W_per_m"] is not None]
    # "realistic" cutoff: keep peak power under ~2 kW/m (a few hundred W for
    # a practical <=0.5 m pipe section -- easily within a resistive
    # cartridge heater or small induction coil's range).
    REALISTIC_POWER_CAP = 2000.0
    candidates = [r for r in feasible if r["P_min_W_per_m"] <= REALISTIC_POWER_CAP]
    chosen = max(candidates, key=lambda r: r["gap_mm"])
    print(f"\nRecommended design: gap={chosen['gap_mm']} mm, frac={chosen['frac']:.2f}, "
          f"P_min={chosen['P_min_W_per_m']:.1f} W/m -> t999={chosen['t999']:.1f} s")

    gap_mm = chosen["gap_mm"]
    frac = chosen["frac"]
    P_design = round(chosen["P_min_W_per_m"] * 1.05, -1)  # +5% margin, rounded

    P = Params(r_wall=0.0127 + gap_mm / 1000.0, n_cells=400, dt=0.1, t_max=400, t_sleeve=0.001)
    grid, is_cu, band = build_case(P, center_frac=frac)
    res = run_full_melt(P, grid, is_cu, P_design, record_every=20,
                         store_snapshots=True, T_setpoint=T_SAFETY_CAP)
    m = res["milestones"]
    print(f"Final verified run (finer grid, +5% power margin = {P_design:.0f} W/m): "
          f"t50={m['t50']:.1f}s t90={m['t90']:.1f}s t99={m['t99']:.1f}s t999={m['t999']:.1f}s")

    # ---- No-sleeve check at the SAME gap, for contrast ----------------------
    grid_bare = Grid(P.r_pipe, P.r_wall, P.n_cells)
    is_cu_bare = np.zeros(grid_bare.n, dtype=bool)
    res_bare = run_full_melt(P, grid_bare, is_cu_bare, 0.0, record_every=20)

    # ---- Save CSVs -----------------------------------------------------------
    with open(f"{RESULTS}/no_sleeve_sweep.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gap_mm", "t999"])
        w.writeheader()
        w.writerows(rows_ns)
    with open(f"{RESULTS}/sleeve_sweep.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gap_mm", "frac", "P_min_W_per_m", "t999"])
        w.writeheader()
        w.writerows(rows_s)

    summary = dict(
        target_s=TARGET_S,
        recommended=dict(
            r_pipe_mm=P.r_pipe * 1000, gap_mm=gap_mm, r_wall_mm=P.r_wall * 1000,
            sleeve_thickness_mm=P.t_sleeve * 1000, sleeve_center_frac=frac,
            sleeve_band_mm=[band[0] * 1000, band[1] * 1000],
            P_design_W_per_m=P_design, T_safety_cap_C=T_SAFETY_CAP,
            milestones=m, no_sleeve_milestones=res_bare["milestones"],
        ),
        no_sleeve_gap_sweep=rows_ns,
        sleeve_gap_sweep=rows_s,
    )
    with open(f"{RESULTS}/summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)

    # ---- Plots --------------------------------------------------------------
    make_plots(P, res, res_bare, rows_ns, rows_s, chosen, band)

    print(f"\nWrote plots, CSVs and summary.json to {RESULTS}/")
    print("\n=== FINAL RECOMMENDATION ===")
    print(f"Pipe radius (fixed, given): {P.r_pipe*1000:.2f} mm")
    print(f"Container inner radius: {P.r_wall*1000:.2f} mm  "
          f"(PCM annulus thickness: {gap_mm} mm, was 32.3 mm)")
    print(f"Copper sleeve: {P.t_sleeve*1000:.1f} mm thick, centered at "
          f"{frac*100:.0f}% of the way from pipe to wall "
          f"({band[0]*1000:.2f}-{band[1]*1000:.2f} mm)")
    print(f"Sleeve pulse power: {P_design:.0f} W/m of pipe length "
          f"(e.g. {P_design*0.3:.0f} W for a 0.3 m section), capped off above "
          f"{T_SAFETY_CAP:.0f} C")
    print(f"Result: fully melted (99.9%) in {m['t999']:.0f} s "
          f"= {m['t999']/60:.2f} min  (target: < 180 s)")
    print(f"Without the sleeve at this same gap: "
          f"{res_bare['milestones']}")


def make_plots(P, res, res_bare, rows_ns, rows_s, chosen, band):
    os.makedirs(RESULTS, exist_ok=True)

    # 1) gap-width sweep: no-sleeve vs with-sleeve required power
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    gaps_ns = [r["gap_mm"] for r in rows_ns]
    t_ns = [r["t999"] if r["t999"] else np.nan for r in rows_ns]
    ax1.plot(gaps_ns, t_ns, "o-", color="steelblue")
    ax1.axhline(180, color="red", linestyle="--", linewidth=1, label="3 min target")
    ax1.set_xlabel("PCM annulus gap width [mm]")
    ax1.set_ylabel("full-melt time [s]")
    ax1.set_title("Pure pipe conduction (no sleeve)")
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.3)

    gaps_s = [r["gap_mm"] for r in rows_s]
    p_s = [r["P_min_W_per_m"] if r["P_min_W_per_m"] else np.nan for r in rows_s]
    ax2.plot(gaps_s, p_s, "o-", color="darkorange")
    ax2.axvline(chosen["gap_mm"], color="green", linestyle=":", label="recommended")
    ax2.set_xlabel("PCM annulus gap width [mm]")
    ax2.set_ylabel("min. sleeve power to melt in <170 s [W/m]")
    ax2.set_title("With pulsed copper sleeve")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{RESULTS}/fig1_gap_sweeps.png", dpi=150)
    plt.close(fig)

    # 2) liquid fraction: recommended design vs no-sleeve at same gap
    fig, ax = plt.subplots(figsize=(7, 4.5))
    h = res["hist"]
    ax.plot(h["t"], h["favg"] * 100, label="with pulsed sleeve", linewidth=1.8)
    hb = res_bare["hist"]
    ax.plot(hb["t"], hb["favg"] * 100, "--", label="no sleeve (same gap)", linewidth=1.8)
    ax.axvline(180, color="red", linestyle=":", linewidth=1, label="3 min target")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("PCM average liquid fraction [%]")
    ax.set_title(f"Melting progress at the recommended {chosen['gap_mm']} mm gap")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{RESULTS}/fig2_liquid_fraction.png", dpi=150)
    plt.close(fig)

    # 3) temperature profile snapshots
    fig, ax = plt.subplots(figsize=(7, 4.5))
    grid, is_cu, _ = build_case(P, center_frac=chosen["frac"])
    snaps = res["snapshots"]
    pick = np.linspace(0, len(snaps) - 1, min(7, len(snaps))).astype(int)
    cmap = plt.cm.viridis(np.linspace(0, 1, len(pick)))
    for c, i in zip(cmap, pick):
        s = snaps[i]
        ax.plot(grid.r_c * 1000, s["T"], color=c, linewidth=1.4, label=f"t={s['t']:.0f} s")
    ax.axhline(P.Tm, color="red", linestyle="--", linewidth=1, label="PCM melt point")
    ax.axvspan(band[0] * 1000, band[1] * 1000, color="orange", alpha=0.25, label="copper sleeve")
    ax.set_xlabel("radius [mm]")
    ax.set_ylabel("temperature [°C]")
    ax.set_title("Radial temperature profiles (recommended fast-melt design)")
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{RESULTS}/fig3_temperature_profiles.png", dpi=150)
    plt.close(fig)

    # 4) sleeve power/temperature trace
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7, 5.5), sharex=True)
    a1.plot(h["t"], h["T_sleeve"], color="darkorange", linewidth=1.4)
    a1.axhline(P.Tm, color="red", linestyle="--", linewidth=0.8, label="PCM melt point")
    a1.axhline(res["hist"]["T_sleeve"].max(), color="gray", linestyle=":", linewidth=0.6)
    a1.set_ylabel("sleeve avg. temp [°C]")
    a1.legend(fontsize=7)
    a1.grid(alpha=0.3)
    a1.set_title("Recommended design: sleeve pulse-heating trace")
    a2.plot(h["t"], h["power"], color="steelblue", linewidth=1.2, drawstyle="steps-post")
    a2.set_xlabel("time [s]")
    a2.set_ylabel("sleeve power [W/m]")
    a2.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{RESULTS}/fig4_sleeve_pulsing.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
