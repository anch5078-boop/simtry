"""
radial_enthalpy_model.py -- Model B: 1-D transient radial enthalpy-method
solver for the PCM annulus surrounding a pulsed vertical heater sleeve,
representative of a single axial cross-section through the heated
section (fixed positions, no sinking -- contrast with Model A).

Gives what Model A cannot: an actual T(r,t) field and a melt-front
r_melt(t) trajectory, including the visible effect of pulsing (melting
stalls -- and the field relaxes -- during each OFF interval).

NUMERICAL METHOD (mirrors PCM_MATLAB Version-1, adapted to 1-D
cylindrical coordinates)
--------------------------------------------------------------------
Finite-volume, backward-Euler in time, harmonic-mean face
conductivities, a few Picard sub-iterations per step for the k(T)/
Capp(T) nonlinearity, and the same *secant* (energy-conserving)
apparent heat capacity used in ``PCM_MATLAB/apparentCapacity.m`` --
that file's docstring explains why the naive tangent-derivative
apparent-capacity form silently under-counts latent heat when a time
step's dT is comparable to the mushy interval width, and how the
secant form fixes it exactly. The same liquid-fraction smoothing
(cubic Hermite / "smoothstep" over [Ts,Tl]) as
``PCM_MATLAB/liquidFraction.m`` is used. The domain is initialized at
P.T_init (well below Ts), not at Tm, for the same reason documented in
PCM_MATLAB/README.md (initializing exactly at Tm falsely puts every
cell at fl=0.5 before it has absorbed any real heat).

Inner boundary at r=r_h: Dirichlet T=T_h while the pulse is ON, and
adiabatic (no term added -- true zero flux, not just "no change") while
OFF. Outer boundary at r=R: always adiabatic (insulated container
wall), matching the Version-1 model's outer-boundary convention.

EFFECTIVE-CONDUCTIVITY SUBMODEL FOR THE MELT LAYER
------------------------------------------------------
Resolving the sub-mm slip-lubricated film and the natural-convection
circulation cell directly on this 1-D grid is out of scope (that is a
multi-dimensional flow problem); instead, once PCM at a location has
melted, its *conductivity* is scaled up by an engineering enhancement
factor E(gap) intended to reproduce the right order of magnitude and
the right qualitative regime crossover, where gap = the current total
melt-layer thickness (r_melt(t) - r_h), not each cell's own depth --
this mirrors how natural-convection "effective conductivity"
correlations (e.g. Raithby & Hollands 1975 for concentric-cylinder
annuli) are normally used: one k_eff for the whole gap, replacing the
conduction profile, not a locally-varying value within it.

    E(gap) = clip( max(slip_term(gap), conv_term(gap)), 1, E_CAP )

    slip_term(gap) = 1 + 3*(b/gap)*exp(-gap/delta_c)
        -- decaying continuation of the flow-rate enhancement factor
        derived in ccm_slip_model.py's Model A (1+3b/delta for
        b<<delta) into the regime where the true dynamic-force-balance
        film no longer applies (delta_c ~ a few hundred um sets where
        that decay happens). This is a deliberate *approximation*: it
        borrows Model A's flow-enhancement algebra as a proxy for
        near-wall convective augmentation of heat transfer, which is a
        reasonable qualitative stand-in but -- unlike Model A's force
        balance -- is not itself separately derived from first
        principles here.

    conv_term(gap) = max(1, 0.386*(Pr/(0.861+Pr))^(1/4) * Ra_gap^(1/4))
        -- standard Raithby-Hollands-type concentric-annulus natural
        convection correlation, Ra_gap = g*beta_l*dT*gap^3/(nu_l*alpha_l),
        dT = T_h - Tm (the characteristic driving temperature
        difference across the melt layer).

    E_CAP = 30 -- both terms are capped; near gap -> 0 the algebraic
        forms are not meant to be extrapolated indefinitely (that
        regime is where Model A's dynamic force balance actually
        governs). 30x is in the range of reported CCM/slip enhancement
        factors over plain conduction and is applied only as a ceiling,
        rarely the binding value away from gap -> 0.

This is explicitly a *simplified engineering blend*, not a
spatially-resolved treatment of two distinct flow mechanisms -- flagged
here rather than left implicit.
"""

from dataclasses import dataclass
import numpy as np

from parameters import Params


def liquid_fraction(T, Ts, Tl):
    """Cubic-Hermite smoothed liquid fraction and its derivative --
    identical formula to PCM_MATLAB/liquidFraction.m."""
    dT = Tl - Ts
    xi = np.clip((T - Ts) / dT, 0.0, 1.0)
    fl = 3 * xi ** 2 - 2 * xi ** 3
    dfldT = 6 * xi * (1 - xi) / dT
    outside = (T <= Ts) | (T >= Tl)
    dfldT = np.where(outside, 0.0, dfldT)
    return fl, dfldT


def apparent_capacity(Tguess, Tn, Ts, Tl, rhocp, rhoL):
    """Secant (energy-conserving) apparent heat capacity -- identical
    method to PCM_MATLAB/apparentCapacity.m; see that file and the
    module docstring above for why this form (rather than the naive
    tangent dfl/dT) is required for energy conservation."""
    dT = Tguess - Tn
    tiny = np.abs(dT) < 1e-8
    flG, _ = liquid_fraction(Tguess, Ts, Tl)
    flN, _ = liquid_fraction(Tn, Ts, Tl)
    with np.errstate(invalid="ignore", divide="ignore"):
        secant = rhocp + rhoL * (flG - flN) / dT
    _, dfldT_tiny = liquid_fraction(Tguess, Ts, Tl)
    Capp = np.where(tiny, rhocp + rhoL * dfldT_tiny, secant)
    return Capp


def enhancement_factor(gap: float, b: float, dT_char: float, P: Params, E_CAP: float = 30.0) -> float:
    """E(gap) as documented in the module docstring."""
    if gap <= 1e-9 or dT_char <= 0:
        return 1.0
    slip_term = 1.0 + 3.0 * (b / gap) * np.exp(-gap / P.delta_c)
    Ra = P.g * P.beta_l * dT_char * gap ** 3 / (P.nu_l * P.alpha_l)
    conv_term = 1.0
    if Ra > 1.0e3:
        Pr = P.Pr
        conv_term = 0.386 * (Pr / (0.861 + Pr)) ** 0.25 * Ra ** 0.25
    E = max(slip_term, conv_term, 1.0)
    return min(E, E_CAP)


def thomas_solve(a, d, c, rhs):
    """Simple Thomas algorithm solving, for each row i,

        a[i]*x[i-1] + d[i]*x[i] + c[i]*x[i+1] = rhs[i]

    i.e. a, c are the *signed* matrix entries (a[0] and c[-1] unused).
    Conductance-form FV assembly (diagonal = +sum of conductances,
    off-diagonals = -conductance to that neighbour) must pass -Gc, not
    +Gc, for a and c -- see the call site in simulate() for the
    conservation argument this sign flip is required for."""
    n = len(d)
    cp = np.empty(n)
    dp = np.empty(n)
    cp[0] = c[0] / d[0]
    dp[0] = rhs[0] / d[0]
    for i in range(1, n):
        m = d[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / m if i < n - 1 else 0.0
        dp[i] = (rhs[i] - a[i] * dp[i - 1]) / m
    x = np.empty(n)
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


@dataclass
class ModelBResult:
    t: np.ndarray
    r_c: np.ndarray             # cell-center radii [m]
    T_history: np.ndarray       # (n_saved, n_r) snapshots
    t_saved: np.ndarray
    r_melt: np.ndarray          # melt-front radius vs time [m]
    liquid_frac_vol: np.ndarray  # volume-weighted overall liquid fraction vs time
    E_used: np.ndarray          # enhancement factor actually applied, vs time
    pulse_on: np.ndarray
    milestones: dict
    b: float


def simulate(P: Params, b: float, save_every: int = 30, verbose: bool = False) -> ModelBResult:
    N = P.n_r
    dr = (P.R - P.r_h) / N
    r_face = P.r_h + np.arange(N + 1) * dr
    r_c = P.r_h + (np.arange(N) + 0.5) * dr
    V = np.pi * (r_face[1:] ** 2 - r_face[:-1] ** 2)  # per unit height

    T = np.full(N, P.T_init)
    dt = P.dt_B
    n_steps = int(P.t_max / dt) + 1

    t_list, r_melt_list, liqfrac_list, E_list, pulse_list = [], [], [], [], []
    T_hist, t_saved = [], []

    dT_char = P.T_h - P.Tm

    picard_iters = 4
    picard_tol = 1e-3

    for step in range(n_steps):
        t_old = step * dt  # computed fresh each step (not accumulated) to
                            # avoid float drift landing on the wrong side of
                            # a pulse edge
        t_new = t_old + dt  # the actual time of the state this step solves
                             # for -- everything recorded below (T, r_melt,
                             # liquid fraction, ...) describes the field at
                             # t_new, not t_old, and must be labeled as such
        is_on = P.pulse_state(t_old)  # BC applied while advancing t_old -> t_new

        # Current melt-front / gap estimate (from the previous step's
        # field) drives this step's enhancement factor; re-evaluated
        # once more after the implicit solve below (one lag is fine --
        # gap changes slowly compared to dt_B).
        fl_prev, _ = liquid_fraction(T, P.Ts, P.Tl)
        gap = dr * np.sum(fl_prev)
        E = enhancement_factor(gap, b, dT_char, P) if is_on else 1.0

        Tn = T.copy()
        Tguess = T.copy()
        for _ in range(picard_iters):
            fl, _ = liquid_fraction(Tguess, P.Ts, P.Tl)
            k = (1.0 - fl) * P.k_s + fl * (P.k_l * E)
            Capp = apparent_capacity(Tguess, Tn, P.Ts, P.Tl, P.rho_s * P.cp, P.rho_s * P.Lf)

            k_face = 2.0 * k[:-1] * k[1:] / (k[:-1] + k[1:])           # interior faces, N-1 of them
            Gc_face = k_face * 2.0 * np.pi * r_face[1:-1] / dr          # interior conductances

            a = np.zeros(N)   # sub-diagonal (signed: -conductance to the i-1 neighbour)
            dgo = np.zeros(N)  # diagonal
            c = np.zeros(N)   # super-diagonal (signed: -conductance to the i+1 neighbour)
            rhs = np.zeros(N)

            cap_term = Capp * V / dt
            dgo[:] = cap_term
            rhs[:] = cap_term * Tn

            # interior face couplings. Row i reads
            #   dgo[i]*T[i] + a[i]*T[i-1] + c[i]*T[i+1] = rhs[i]
            # with dgo[i] = cap_term[i] + (sum of neighbouring Gc_face) and
            # a[i]/c[i] = -Gc_face (moving the neighbour term to the LHS
            # flips its sign) -- passing +Gc_face here instead silently
            # breaks energy conservation (thomas_solve then solves a
            # different, non-conservative system) even though it still
            # looks like a converged, well-behaved solve.
            a[1:] -= Gc_face
            c[:-1] -= Gc_face
            dgo[:-1] += Gc_face
            dgo[1:] += Gc_face

            # inner wall (r=r_h): Dirichlet T_h while ON, adiabatic while OFF
            if is_on:
                Gc_wall = k[0] * 2.0 * np.pi * P.r_h / (0.5 * dr)
                dgo[0] += Gc_wall
                rhs[0] += Gc_wall * P.T_h
            # outer wall (r=R): always adiabatic -- no term added

            Tnew = thomas_solve(a, dgo, c, rhs)
            change = np.max(np.abs(Tnew - Tguess))
            scale = max(1.0, np.max(np.abs(Tguess)))
            Tguess = Tnew
            if change < picard_tol * scale:
                break

        T = Tguess
        fl_now, _ = liquid_fraction(T, P.Ts, P.Tl)

        # melt front: linear-interpolated crossing of T=Tm, same
        # convention as PCM_MATLAB/frontPosition.m
        above = np.where(T > P.Tm)[0]
        if above.size == 0:
            r_melt = r_c[0]
        elif above[-1] == N - 1:
            r_melt = r_c[-1]
        else:
            idx = above[-1]
            T1, T2 = T[idx], T[idx + 1]
            r1, r2 = r_c[idx], r_c[idx + 1]
            frac = (P.Tm - T1) / (T2 - T1)
            r_melt = r1 + frac * (r2 - r1)

        liquid_frac_vol = float(np.sum(fl_now * V) / np.sum(V))

        t_list.append(t_new)
        r_melt_list.append(r_melt)
        liqfrac_list.append(liquid_frac_vol)
        E_list.append(E)
        pulse_list.append(is_on)

        if step % save_every == 0:
            T_hist.append(T.copy())
            t_saved.append(t_new)

        if liquid_frac_vol > 0.999:
            break

    t_arr = np.array(t_list)
    liqfrac_arr = np.array(liqfrac_list)
    milestones = {}
    for label, thresh in (("t50", 0.50), ("t90", 0.90), ("t95", 0.95), ("t99", 0.99)):
        idx = np.argmax(liqfrac_arr >= thresh) if np.any(liqfrac_arr >= thresh) else -1
        if idx <= 0:
            milestones[label] = None
        else:
            f0, f1 = liqfrac_arr[idx - 1], liqfrac_arr[idx]
            t0, t1 = t_arr[idx - 1], t_arr[idx]
            frac = 0.0 if f1 == f0 else (thresh - f0) / (f1 - f0)
            milestones[label] = t0 + frac * (t1 - t0)

    if verbose:
        print(f"[Model B, b={b*1e6:.0f} um] t50={milestones['t50']}, "
              f"t90={milestones['t90']}, t99={milestones['t99']}")

    return ModelBResult(
        t=t_arr, r_c=r_c,
        T_history=np.array(T_hist), t_saved=np.array(t_saved),
        r_melt=np.array(r_melt_list), liquid_frac_vol=liqfrac_arr,
        E_used=np.array(E_list), pulse_on=np.array(pulse_list),
        milestones=milestones, b=b,
    )
