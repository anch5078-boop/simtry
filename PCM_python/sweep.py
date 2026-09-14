"""
sweep.py -- vary the copper sleeve's radial position across the
pipe->wall gap (holding the "lowest power to reach 40 C within the
ramp target" control rule fixed at each position) to find the
placement that minimizes total melting time.
"""

from __future__ import annotations

import numpy as np

from model import Grid
from params import Params
from simulate import build_case, find_min_power, run_full_melt


def run_position_sweep(P: Params, fracs=None, record_every=60, verbose=True):
    if fracs is None:
        fracs = np.linspace(0.15, 0.85, 8)

    rows = []
    for frac in fracs:
        grid, is_cu, (r_in, r_out) = build_case(P, center_frac=frac)
        P_min, T_end = find_min_power(P, grid, is_cu)
        result = run_full_melt(P, grid, is_cu, P_min, record_every=record_every)
        m = result["milestones"]
        row = dict(
            frac=frac, r_in=r_in, r_out=r_out,
            r_mid=0.5 * (r_in + r_out),
            P_min_W_per_m=P_min,
            t50=m["t50"], t90=m["t90"], t99=m["t99"], t999=m["t999"],
            duty_cycle=result["duty_cycle"],
            converged=result["converged"],
        )
        rows.append(row)
        if verbose:
            def fmt(x):
                return f"{x:8.1f}" if x is not None else "    n/a "
            print(f"  frac={frac:4.2f}  r_mid={row['r_mid']*1000:6.2f}mm  "
                  f"P_min={P_min:7.2f} W/m  "
                  f"t50={fmt(m['t50'])}s  t90={fmt(m['t90'])}s  "
                  f"t99={fmt(m['t99'])}s  t999={fmt(m['t999'])}s")
    return rows


def run_no_sleeve_reference(P: Params, record_every=60):
    """Pipe-only reference: no copper sleeve anywhere (pure PCM from
    pipe to wall) -- shows the benefit (or not) of the auxiliary
    sleeve heater."""
    grid = Grid(P.r_pipe, P.r_wall, P.n_cells)
    is_cu = np.zeros(grid.n, dtype=bool)
    result = run_full_melt(P, grid, is_cu, 0.0, record_every=record_every)
    return grid, result
