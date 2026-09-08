"""
ccm_slip_model.py -- Model A: 0-D lumped force-balance / lubrication
model for slip-enhanced close-contact melting (sCCM) of a solid PCM
column sinking under gravity around a vertical, pulsed cylindrical
heater sleeve.

PHYSICAL PICTURE
-----------------
The tube is vertical, so unlike the horizontal grooved-plate CCM setup
in the reference SI (gravity normal to the melting plate), here gravity
acts *along* the tube axis. The heater is a sleeve running up the axis;
solid PCM fills the annulus around it and, over the whole column
height, well above the heated section too. As PCM melts off the
sleeve's surface, a thin lubricating melt film separates the sleeve
from the (still solid, still rigid) PCM column above it. The column's
own net weight (gravity minus buoyancy) presses down on this film and
squeezes melt axially downward out of the heated section; mass
conservation (melting flux in = axial drainage flux out) plus a
lubrication force balance sets the film's quasi-steady thickness
delta(t), from which the heat flux q = k_l*dT/delta follows. This is
the vertical-tube analogue of classic horizontal CCM (Bejan;
Moallemi & Viskanta 1985) and reduces to it algebraically -- see
"Sanity check" below.

The slip coating (Li et al. 2026) does not change the conduction
resistance across the film directly (heat still crosses by conduction,
q = k_l*dT/delta, exactly as without slip) -- it changes the
*equilibrium* delta the force balance settles on, by reducing the
viscous resistance to the axial drainage flow that must be squeezed
out through the film. A thinner equilibrium film means a higher
q. This matches the mechanism reported in the reference paper.

DERIVATION (plane-unrolled lubrication, one slip wall)
-------------------------------------------------------
Unroll the annular film (width delta << r_h) into a 2-D channel of
gap delta, spanwise width w = 2*pi*r_h, running the wetted height Hc.
Plane Poiseuille flow driven by -dp/dz = G, with Navier slip length b
on the heater wall (y=0) and no-slip on the solid's melting surface
(y=delta):

    u(y) = -(G/2mu) y^2 + C1 y + C2,   C1 = G delta^2 / (2 mu (delta+b)),
    C2 = b C1

Integrating u over the gap gives the flow rate per unit width

    Q' = G * delta^3 (delta + 4b) / (12 mu (delta + b))
       = G * delta_eff^3 / (12 mu),   delta_eff^3 := delta^3(delta+4b)/(delta+b)

which recovers the ordinary no-slip result Q' = G delta^3/(12 mu) at
b=0, and the b -> infinity (shear-free wall) limit Q' = G delta^3/(3 mu)
(4x enhancement), both of which are standard checks on a one-wall-slip
plane Poiseuille solution.

Melting is distributed uniformly over the wetted height Hc (uniform-dT,
uniform-delta lumped approximation), so the axial volumetric flow at
height z (measured from the open bottom, z=0, up to the closed top of
the wetted zone, z=Hc) is the melt generated above it:

    Q(z) = q_h w (Hc - z) / (rho_l Lf),   q_h = k_l dT / delta

Using dp/dz = 12 mu Q(z) / (w delta_eff^3) and integrating twice gives
the excess pressure profile p(z), and integrating p(z) over the
contact area gives the total upward force supporting the column's net
weight W_net:

    W_net = 4 mu_l w q_h Hc^3 / (delta_eff^3 rho_l Lf)                  (*)

With q_h = k_l dT / delta substituted, (*) is one nonlinear equation in
the one unknown delta (given W_net, dT, Hc, b); solved numerically
per time step with a bisection root-finder.

Sanity check (b=0):  delta * delta_eff^3 = delta^4, so
    delta = [ 4 mu_l w k_l dT Hc^3 / (W_net rho_l Lf) ]^(1/4)
which is exactly the classic CCM 1/4-power scaling for film thickness
(propto (mu k dT L^3 / (W rho_L))^(1/4)) reported throughout the CCM
literature (Bejan 1994; Moallemi & Viskanta 1985) -- confirms the
lumped derivation above is dimensionally and structurally consistent
with the established (no-slip) theory before slip is added.

MODEL ASSUMPTIONS / LIMITATIONS (stated explicitly, not hidden)
------------------------------------------------------------------
  * Lumped in z: dT and delta are treated as spatially uniform over the
    wetted height Hc at each instant (only the *flow rate* Q(z) is
    allowed to vary with z, which is what the force balance needs).
    A fully-resolved sCCM model would let delta vary with z; this
    lumped version is a standard first simplification in the CCM
    literature and is adequate for overall charge-time trends.
  * Pulsing: the thin film's own thermal mass is neglected, so during
    an OFF interval (wall adiabatic) dT -> 0 instantly and melting
    pauses (q_h=0) rather than decaying gradually. This is conservative
    in the sense that it does not credit the model with any "coasting"
    melting during OFF.
  * The force balance is algebraic (quasi-steady), not an ODE for
    delta itself, so delta(t=0) is simply whatever solves (*) at the
    initial (finite) W_net -- there is no zero-gap start-up
    singularity to regularize. What the model *does* skip over is the
    brief real conduction-limited transient before quasi-steady
    lubrication is established (well documented in the CCM
    literature, e.g. Moallemi & Viskanta 1985); this is short compared
    to the charge times of interest here and is not separately
    modeled. ``delta0`` instead serves as a minimum achievable film
    thickness (surface roughness / non-ideal-contact floor): delta is
    clamped to it, which in turn caps q_h wherever the smooth-wall
    force balance would otherwise want an unrealistically thin film.
  * Breaks down (delta -> unbounded) as W_net -> 0 near complete
    melt-out (Hs -> 0); capped at solve_delta's upper ``delta_bounds``
    and flagged in the returned arrays via ``capped``.
"""

from dataclasses import dataclass
import math
import numpy as np
from scipy.optimize import brentq

from parameters import Params


def delta_eff_cubed(delta: np.ndarray, b: float) -> np.ndarray:
    """delta_eff^3 = delta^3 (delta + 4b) / (delta + b) -- the
    slip-modified effective hydrodynamic gap (cubed) for one-wall
    Navier-slip plane Poiseuille flow. Reduces to delta**3 at b=0."""
    return delta ** 3 * (delta + 4.0 * b) / (delta + b)


def solve_delta(Wnet: float, mu_l: float, k_l: float, dT: float, Hc: float,
                 w: float, rho_l: float, Lf: float, b: float,
                 delta_bounds=(1e-9, 5e-2)) -> float:
    """Root-find the quasi-steady film thickness delta [m] satisfying
    the force balance (*), given the current driving dT and wetted
    height Hc. Returns delta_bounds[1] (capped) if no root exists in
    range (i.e. even the widest allowed film can't be pushed
    thin/thick enough to balance -- occurs only as Wnet -> 0)."""
    lo, hi = delta_bounds

    def predicted_Wnet(delta):
        return 4.0 * mu_l * w * Hc ** 3 * k_l * dT / (delta * delta_eff_cubed(delta, b) * rho_l * Lf)

    f_lo = predicted_Wnet(lo) - Wnet
    f_hi = predicted_Wnet(hi) - Wnet
    if f_lo <= 0:
        # Even the thinnest allowed film can't support this much
        # weight at this dT -- physically means the column would need
        # to be arrested by something else (not modeled); return the
        # thinnest bound as the best available estimate.
        return lo
    if f_hi >= 0:
        # Weight is small enough that even the widest allowed film
        # would "support" it -- near complete melt-out. Cap.
        return hi
    return brentq(lambda d: predicted_Wnet(d) - Wnet, lo, hi, xtol=1e-12, rtol=1e-10)


@dataclass
class ModelAResult:
    t: np.ndarray
    Hs: np.ndarray          # remaining solid column height [m]
    melt_fraction: np.ndarray
    delta: np.ndarray       # film thickness [m] (nan while OFF)
    q_h: np.ndarray         # heater-surface heat flux [W/m^2]
    V_sink: np.ndarray      # solid descent (~melt) velocity [m/s]
    pulse_on: np.ndarray    # bool
    capped: np.ndarray      # bool, True where delta hit a solver bound
    milestones: dict        # {"t50": ..., "t90": ..., "t95": ..., "t99": ...} seconds or None
    b: float


def _milestones_from_series(t: np.ndarray, melt_fraction: np.ndarray) -> dict:
    # t99 is included for parity with the Version-1 MATLAB milestone
    # convention, but note the caveat in this module's docstring: as
    # the remaining solid (and hence its driving weight W_net) shrinks
    # towards zero, delta is pushed towards its upper solver bound and
    # q_h towards zero, so the last percent of melt can have a long,
    # slowly-converging tail (or not converge within t_max at all).
    # t95 is reported alongside it as the more robust summary number
    # for this model.
    out = {}
    for label, thresh in (("t50", 0.50), ("t90", 0.90), ("t95", 0.95), ("t99", 0.99)):
        idx = np.argmax(melt_fraction >= thresh) if np.any(melt_fraction >= thresh) else -1
        if idx <= 0:
            out[label] = None
            continue
        # linear interpolation between idx-1 and idx for a smoother estimate
        f0, f1 = melt_fraction[idx - 1], melt_fraction[idx]
        t0, t1 = t[idx - 1], t[idx]
        frac = 0.0 if f1 == f0 else (thresh - f0) / (f1 - f0)
        out[label] = t0 + frac * (t1 - t0)
    return out


def simulate(P: Params, b: float, verbose: bool = False) -> ModelAResult:
    """Explicit time march of the lumped sCCM ODE for one slip length b.

    State variable: Hs, the height of solid PCM column remaining above
    (and including, down to) the melting interface. Hs(0) = P.H_col.
    """
    w = 2.0 * math.pi * P.r_h
    A_cross = P.A_cross
    dt = P.dt_A
    n_steps = int(P.t_max / dt) + 1
    Hs_floor = 1.0e-4  # treat as "fully melted" below this remaining height

    t_arr = np.zeros(n_steps)
    Hs_arr = np.zeros(n_steps)
    delta_arr = np.full(n_steps, np.nan)
    qh_arr = np.zeros(n_steps)
    Vsink_arr = np.zeros(n_steps)
    pulse_arr = np.zeros(n_steps, dtype=bool)
    capped_arr = np.zeros(n_steps, dtype=bool)

    Hs = P.H_col
    n_used = n_steps
    for i in range(n_steps):
        t = i * dt  # computed fresh each step (not accumulated) to avoid
                    # float drift landing on the wrong side of a pulse edge
        t_arr[i] = t
        Hs_arr[i] = Hs

        is_on = P.pulse_state(t)
        pulse_arr[i] = is_on
        Hc = min(Hs, P.H_heater)

        if is_on and Hc > 0.0 and Hs > Hs_floor:
            dT = P.T_h - P.Tm
            Wnet = (P.rho_s - P.rho_l) * P.g * A_cross * Hs
            if Wnet > 0.0:
                delta = solve_delta(Wnet, P.mu_l, P.k_l, dT, Hc, w, P.rho_l, P.Lf, b)
                capped_arr[i] = (delta <= 1e-9 * 1.0001) or (delta >= 5e-2 * 0.9999)
                delta = max(delta, P.delta0)  # surface-roughness / non-ideal-contact floor
                q_h = P.k_l * dT / delta
            else:
                delta, q_h = np.nan, 0.0
        else:
            delta, q_h = np.nan, 0.0

        delta_arr[i] = delta
        qh_arr[i] = q_h

        V_sink = q_h * (2.0 * math.pi * P.r_h * Hc) / (P.rho_s * P.Lf * A_cross) if Hc > 0 else 0.0
        Vsink_arr[i] = V_sink

        if Hs <= Hs_floor:
            n_used = i + 1
            break

        Hs = max(Hs - V_sink * dt, 0.0)

    t_arr = t_arr[:n_used]
    Hs_arr = Hs_arr[:n_used]
    delta_arr = delta_arr[:n_used]
    qh_arr = qh_arr[:n_used]
    Vsink_arr = Vsink_arr[:n_used]
    pulse_arr = pulse_arr[:n_used]
    capped_arr = capped_arr[:n_used]

    melt_fraction = 1.0 - Hs_arr / P.H_col
    milestones = _milestones_from_series(t_arr, melt_fraction)

    if verbose:
        print(f"[Model A, b={b*1e6:.0f} um] t50={milestones['t50']}, "
              f"t90={milestones['t90']}, t99={milestones['t99']}")

    return ModelAResult(t=t_arr, Hs=Hs_arr, melt_fraction=melt_fraction,
                         delta=delta_arr, q_h=qh_arr, V_sink=Vsink_arr,
                         pulse_on=pulse_arr, capped=capped_arr,
                         milestones=milestones, b=b)
