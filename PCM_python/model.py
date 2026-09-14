"""
model.py -- 1-D radial (cylindrical) transient conduction + phase-change
solver for the "hot pipe / PCM / pulsed copper sleeve / PCM / insulated
container wall" geometry.

Numerical method (deliberately mirrors PCM_MATLAB/, the validated
Version-1 solver already in this repo -- see PCM_MATLAB/README.md for
the reasoning; ported here to a 1-D cylindrical grid because this
geometry is axisymmetric and uniform along the pipe axis, so a full 3-D
Cartesian grid is unnecessary):

  - Governing equation, cylindrical form (per unit axial length):
        C_app(T) dT/dt = (1/r) d/dr( r k(T) dT/dr ) + q'''
  - Smoothed cubic-Hermite liquid fraction f_l(T) over a narrow mushy
    interval [Ts, Tl] straddling the nominal melting point Tm.
  - Energy-conserving *secant* apparent heat capacity (not the naive
    tangent dfl/dT form -- see PCM_MATLAB/apparentCapacity.m for why
    that under-counts latent heat when a step's dT is comparable to the
    mushy width). This is the single most important correctness
    detail in the whole scheme.
  - Backward-Euler implicit time stepping with a few Picard
    sub-iterations per step to handle the k(T)/C_app(T) nonlinearity.
  - Harmonic-mean face conductivity (correct for two materials in
    series, e.g. the PCM/copper interface where k differs ~2000x).
  - Inner boundary: Dirichlet (fixed T = hot-water temperature) at the
    pipe outer surface -- i.e. the pipe wall itself is assumed to have
    negligible thermal resistance compared to the PCM annulus, so the
    PCM-facing surface simply sits at the hot-water temperature. This
    is the standard simplifying assumption for a first screening model
    (thin metal pipe wall, well-mixed turbulent water -> high internal
    h); note it in any write-up as a modelling choice.
  - Outer boundary: adiabatic (perfectly insulated container wall) --
    there is no other heat sink in this problem, matching the
    "container" in the diagram.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import solve_banded


# --------------------------------------------------------------------------
# Smoothed liquid fraction (cubic Hermite / smoothstep) -- same functional
# form as PCM_MATLAB/liquidFraction.m
# --------------------------------------------------------------------------
def liquid_fraction(T, Ts, Tl):
    dT = Tl - Ts
    xi = (T - Ts) / dT
    xi = np.clip(xi, 0.0, 1.0)
    fl = 3 * xi**2 - 2 * xi**3
    return fl


def liquid_fraction_deriv(T, Ts, Tl):
    dT = Tl - Ts
    xi = np.clip((T - Ts) / dT, 0.0, 1.0)
    dfldT = 6 * xi * (1 - xi) / dT
    outside = (T <= Ts) | (T >= Tl)
    dfldT = np.where(outside, 0.0, dfldT)
    return dfldT


def apparent_capacity(Tguess, Tn, Ts, Tl, rhocp, rhoL):
    """Energy-conserving secant apparent heat capacity -- see
    PCM_MATLAB/apparentCapacity.m for the derivation. Falls back to the
    tangent derivative where Tguess==Tn to machine precision."""
    dT = Tguess - Tn
    tiny = np.abs(dT) < 1e-8
    Capp = np.empty_like(Tguess, dtype=float)

    if np.any(~tiny):
        flG = liquid_fraction(Tguess[~tiny], Ts, Tl)
        flN = liquid_fraction(Tn[~tiny], Ts, Tl)
        Capp[~tiny] = rhocp + rhoL * (flG - flN) / dT[~tiny]

    if np.any(tiny):
        dfldT = liquid_fraction_deriv(Tguess[tiny], Ts, Tl)
        Capp[tiny] = rhocp + rhoL * dfldT

    return Capp


# --------------------------------------------------------------------------
# Geometry / grid
# --------------------------------------------------------------------------
class Grid:
    """Uniform-spacing 1-D radial finite-volume grid from r_in to r_out.
    Per-unit-axial-length quantities throughout (so 'volume' is really an
    area [m^2] and 'face area' is really a circumference [m]) -- this is
    the standard reduction of an axisymmetric, axially-uniform problem to
    one radial dimension."""

    def __init__(self, r_in, r_out, n_cells):
        self.n = n_cells
        self.r_face = np.linspace(r_in, r_out, n_cells + 1)
        self.r_c = 0.5 * (self.r_face[:-1] + self.r_face[1:])
        self.dr = self.r_face[1] - self.r_face[0]
        self.V = np.pi * (self.r_face[1:] ** 2 - self.r_face[:-1] ** 2)  # m^2 (per unit length)
        self.A_face = 2 * np.pi * self.r_face  # m (per unit length), length n+1


def make_material_mask(grid: Grid, r_sleeve_in, r_sleeve_out):
    """True where the cell centre falls inside the copper sleeve band."""
    return (grid.r_c >= r_sleeve_in) & (grid.r_c <= r_sleeve_out)


# --------------------------------------------------------------------------
# One implicit time step (backward Euler + Picard), Dirichlet inner BC,
# adiabatic outer BC.
# --------------------------------------------------------------------------
def solve_one_step(Tn, grid: Grid, is_cu, P_sleeve_watts_per_m, dt, mat, picard_iters=4, picard_tol=1e-3):
    """
    Tn                : (N,) temperature at start of step [C]
    is_cu             : (N,) bool mask, True in copper-sleeve cells
    P_sleeve_watts_per_m : total instantaneous sleeve power this step [W/m]
    mat               : dict of material properties (see params.py)
    """
    N = grid.n
    V = grid.V
    dr = grid.dr
    is_pcm = ~is_cu

    V_cu = V[is_cu].sum()
    q_cell = np.zeros(N)
    if V_cu > 0 and P_sleeve_watts_per_m != 0.0:
        q_cell[is_cu] = P_sleeve_watts_per_m * V[is_cu] / V_cu  # distribute by volume

    Tguess = Tn.copy()
    Tnew = Tn.copy()

    for _ in range(picard_iters):
        fl = liquid_fraction(Tguess, mat["Ts"], mat["Tl"])
        k = np.where(is_pcm, (1 - fl) * mat["k_pcm_s"] + fl * mat["k_pcm_l"], mat["k_cu"])

        Capp = np.empty(N)
        Capp[is_pcm] = apparent_capacity(
            Tguess[is_pcm], Tn[is_pcm], mat["Ts"], mat["Tl"], mat["rhocp_pcm"], mat["rhoL_pcm"]
        )
        Capp[is_cu] = mat["rhocp_cu"]

        # Harmonic-mean interior face conductance (uniform dr -> centres
        # equidistant from the shared face).
        k_face = 2 * k[:-1] * k[1:] / (k[:-1] + k[1:])
        G_int = k_face * grid.A_face[1:-1] / dr  # length N-1

        # Inner Dirichlet boundary: pipe surface at r_face[0], distance
        # dr/2 to the first cell centre.
        G_bc_in = 2 * k[0] * grid.A_face[0] / dr

        diag = Capp * V / dt
        diag[:-1] += G_int
        diag[1:] += G_int
        diag[0] += G_bc_in

        off = -G_int

        # Tridiagonal system in LAPACK banded storage (solve_banded is far
        # cheaper here than a general sparse solve -- this runs many
        # thousands of times per case, and many cases per sweep).
        ab = np.zeros((3, N))
        ab[0, 1:] = off      # super-diagonal
        ab[1, :] = diag      # main diagonal
        ab[2, :-1] = off     # sub-diagonal

        b = Capp * V / dt * Tn + q_cell
        b[0] += G_bc_in * mat["T_hotwater"]

        Tnew = solve_banded((1, 1), ab, b)

        change = np.max(np.abs(Tnew - Tguess))
        scale = max(1.0, np.max(np.abs(Tguess)))
        Tguess = Tnew
        if change < picard_tol * scale:
            break

    fl_final = liquid_fraction(Tnew, mat["Ts"], mat["Tl"])
    return Tnew, fl_final
