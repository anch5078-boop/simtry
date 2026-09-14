"""
simulate.py -- build a case (grid + sleeve position), find the minimum
sleeve peak power that reaches the PCM melt point within the target
ramp time, then run the full transient melt simulation with a bang-bang
thermostat holding the sleeve at the melt point thereafter.
"""

from __future__ import annotations

import numpy as np

from model import Grid, make_material_mask, solve_one_step, liquid_fraction
from params import Params


def build_case(P: Params, center_frac=None):
    grid = Grid(P.r_pipe, P.r_wall, P.n_cells)
    r_in, r_out = P.sleeve_radii(center_frac)
    is_cu = make_material_mask(grid, r_in, r_out)
    if not np.any(is_cu):
        raise ValueError(
            f"No grid cell falls inside the sleeve band [{r_in}, {r_out}] m -- "
            f"increase n_cells or t_sleeve."
        )
    return grid, is_cu, (r_in, r_out)


def _volume_avg(field, V, mask):
    return float(np.sum(field[mask] * V[mask]) / np.sum(V[mask]))


def ramp_final_sleeve_temp(P: Params, grid: Grid, is_cu, P_peak, t_ramp):
    """Run the sleeve-heating ramp (heater full-on at P_peak from a cold
    start) for t_ramp seconds and return the volume-averaged sleeve
    temperature at the end. Used as the scalar objective for the
    minimum-power bisection."""
    mat = P.material_dict()
    N = grid.n
    T = np.full(N, P.T_init)
    n_steps = max(1, int(round(t_ramp / P.dt)))
    for _ in range(n_steps):
        T, _ = solve_one_step(T, grid, is_cu, P_peak, P.dt, mat,
                               P.picard_iters, P.picard_tol)
    return _volume_avg(T, grid.V, is_cu)


def find_min_power(P: Params, grid: Grid, is_cu, p_lo=0.0, p_hi=None,
                    tol_temp=0.05, max_iter=40):
    """Bisect for the smallest sleeve peak power [W/m of pipe length]
    such that the volume-averaged sleeve temperature reaches P.Tm within
    P.t_ramp_target seconds of a cold start (with the pipe boundary
    already at T_hotwater throughout, as in the full problem)."""
    if p_hi is None:
        # Grow p_hi until it clearly overshoots Tm within the ramp window.
        p_hi = 50.0
        for _ in range(30):
            T_end = ramp_final_sleeve_temp(P, grid, is_cu, p_hi, P.t_ramp_target)
            if T_end >= P.Tm:
                break
            p_hi *= 2.0
        else:
            raise RuntimeError("Could not bracket a sufficient sleeve power.")

    T_lo = ramp_final_sleeve_temp(P, grid, is_cu, p_lo, P.t_ramp_target)
    if T_lo >= P.Tm:
        return p_lo, T_lo  # even zero power gets there (e.g. very slow ramp target)

    lo, hi = p_lo, p_hi
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        T_end = ramp_final_sleeve_temp(P, grid, is_cu, mid, P.t_ramp_target)
        if T_end >= P.Tm:
            hi = mid
        else:
            lo = mid
        if abs(T_end - P.Tm) < tol_temp:
            break
    return hi, T_end


def run_full_melt(P: Params, grid: Grid, is_cu, P_peak, record_every=30,
                   store_snapshots=False, T_setpoint=None):
    """Full transient simulation: bang-bang thermostat holds the sleeve
    at T_setpoint (heater ON at P_peak while volume-avg sleeve T <
    T_setpoint, OFF otherwise) while the pipe boundary sits fixed at
    T_hotwater throughout. T_setpoint defaults to P.Tm (the efficiency-
    oriented "lowest power to just reach/hold the melt point" rule);
    pass a higher value (or np.inf, i.e. always-on) for a speed-oriented
    "run the sleeve as hard as the supply allows" case. Returns a dict
    of time histories and t50/t90/t99/t999 (99.9%) melt milestones based
    on the volume-weighted average PCM liquid fraction. If
    store_snapshots, also returns full T(r) and fl(r) profiles at each
    recorded step (for temperature-profile plots)."""
    if T_setpoint is None:
        T_setpoint = P.Tm
    mat = P.material_dict()
    N = grid.n
    is_pcm = ~is_cu
    V = grid.V

    has_sleeve = bool(np.any(is_cu))

    T = np.full(N, P.T_init)
    t = 0.0
    n_steps = int(round(P.t_max / P.dt))

    hist = dict(t=[0.0], favg=[0.0], T_sleeve=[P.T_init if has_sleeve else float("nan")],
                T_max=[P.T_init], power=[0.0], duty_on=[False])
    snapshots = [dict(t=0.0, T=T.copy(), fl=np.zeros(N))] if store_snapshots else None

    milestones = {"t50": None, "t90": None, "t99": None, "t999": None}
    heater_on_time = 0.0

    for step in range(1, n_steps + 1):
        if has_sleeve:
            T_sleeve_now = _volume_avg(T, V, is_cu)
            on = T_sleeve_now < T_setpoint
        else:
            on = False
        P_applied = P_peak if on else 0.0
        if on:
            heater_on_time += P.dt

        T, fl = solve_one_step(T, grid, is_cu, P_applied, P.dt, mat,
                                P.picard_iters, P.picard_tol)
        t += P.dt

        favg = _volume_avg(fl, V, is_pcm)

        if step % record_every == 0:
            hist["t"].append(t)
            hist["favg"].append(favg)
            hist["T_sleeve"].append(_volume_avg(T, V, is_cu) if has_sleeve else float("nan"))
            hist["T_max"].append(float(np.max(T)))
            hist["power"].append(P_applied)
            hist["duty_on"].append(on)
            if store_snapshots:
                snapshots.append(dict(t=t, T=T.copy(), fl=fl.copy()))

        for name, thresh in (("t50", 0.50), ("t90", 0.90), ("t99", 0.99), ("t999", 0.999)):
            if milestones[name] is None and favg >= thresh:
                milestones[name] = t

        if milestones["t999"] is not None:
            break

    for k in hist:
        hist[k] = np.array(hist[k])

    duty_cycle_overall = heater_on_time / t if t > 0 else float("nan")

    return dict(hist=hist, milestones=milestones, t_final=t,
                duty_cycle=duty_cycle_overall, converged=milestones["t999"] is not None,
                favg_final=float(hist["favg"][-1]), snapshots=snapshots)
