"""
Close-contact melting (CCM) on structured slip surfaces.

This module implements a slip-modified lubrication model of close-contact
melting, built to capture the same physical mechanism as:

    S. Hu, N. Hu, Y. Lai, Z. Li, X. Gao, L. Fan,
    "Close-contact melting regulated by structured slip surfaces,"
    Appl. Phys. Lett. 128, 073903 (2026). https://doi.org/10.1063/5.0311034

WHAT THIS IS AND ISN'T
-----------------------
The paper's exact reduced-order model (its Eqs. 2-3) is a nondimensional
two-scale system whose closure -- the explicit formulas for the flat-meniscus
slip lengths lambda^(0), lambda_t^(0) and the meniscus-curvature correction
lambda^(1), lambda_t^(1) for a given groove period / gas fraction -- lives in
that paper's Supplementary Material (Secs. S1-S8) and an earlier companion
paper (their Ref. 23), neither of which is available to this implementation
(only the 7-page main article was supplied).

Rather than guess at the exact, unpublished nondimensional PDE coefficients,
this module re-derives a physically equivalent model from first principles:

  1. Macroscopic Stefan (energy) balance at the melting front, augmented by
     an effective *thermal* slip length (the film behaves thermally as if it
     were (h + lambda_t) thick).
  2. A slip-modified 1-D lubrication (Reynolds) equation for the film
     pressure P(x, t) along the drainage direction, which is what produces
     the *spatial* variation of pressure -- and hence of meniscus curvature
     -- that is the paper's central point.
  3. A force balance: the integral of film pressure supports the current
     weight of the unmelted solid, which self-consistently sets the (spatially
     averaged) film thickness at each instant.
  4. A meniscus-curvature feedback on the local effective slip length: locally
     higher film pressure bows the gas-liquid meniscus further into the
     groove (Young-Laplace, R = sigma / P), which *reduces* the local slip
     length relative to its flat-meniscus value -- reproducing the paper's
     key qualitative finding that the meniscus-resolving ("refined") model
     melts *slower* than the constant-slip ("conventional") model.

The flat-meniscus slip length lambda^(0) uses the well-established closed
form for longitudinal micro-ridge / superhydrophobic striped surfaces
(Philip, J. Appl. Math. Phys. 1972; used in this exact form by Lauga & Stone,
J. Fluid Mech. 2003), which is standard, citable, and independent of the
paper's unpublished supplementary formulas. The curvature-feedback closure
(``MENISCUS_COEFF`` below) is a documented, tunable, first-order
approximation -- it is the one piece of this model that is a deliberate
stand-in for the paper's proprietary closure, and is flagged as such
everywhere it is used.
"""

from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq
from scipy.linalg import solve_banded


# ---------------------------------------------------------------------------
# Physical properties: ice / water at ~0 degC (standard reference values,
# consistent with the paper's use of ice as the representative PCM).
# ---------------------------------------------------------------------------
RHO_S = 917.0        # kg/m^3, ice density
RHO_L = 1000.0       # kg/m^3, water density (melt film)
MU_L = 1.5e-3        # Pa s,  water dynamic viscosity (film-averaged, 0-4 degC)
K_L = 0.56           # W/(m K), water thermal conductivity
LATENT_HEAT = 3.34e5  # J/kg,  latent heat of fusion of ice
CP_L = 4186.0        # J/(kg K), water specific heat
SIGMA = 0.0756       # N/m,   water-air surface tension near 0 degC
GRAVITY = 9.81       # m/s^2


@dataclass
class Groove:
    """Longitudinal micro-groove geometry on the heated plate."""
    pitch: float   # l*, groove period [m]
    phi: float     # gas (shear-free) fraction of one period, 0 <= phi < 1
    name: str = ""

    @property
    def lambda0(self) -> float:
        """Flat-meniscus (leading-order) effective slip length [m].

        Philip (1972) / Lauga & Stone (2003) closed form for longitudinal
        micro-ridges in the deep-channel limit:
            lambda0 = (pitch/pi) * ln( sec(pi*phi/2) )
        """
        if self.phi <= 0.0:
            return 0.0
        phi = min(self.phi, 0.999)
        return (self.pitch / np.pi) * np.log(1.0 / np.cos(np.pi * phi / 2))


FLAT = Groove(pitch=1.0, phi=0.0, name="Smooth (no-slip)")


def pitch_for_lambda0(lambda0_target: float, phi: float) -> float:
    """Invert the Philip/Lauga-Stone formula: the groove pitch that gives a
    desired flat-meniscus slip length at a given gas fraction phi."""
    phi = min(phi, 0.999)
    return lambda0_target * np.pi / np.log(1.0 / np.cos(np.pi * phi / 2))


# ---------------------------------------------------------------------------
# Meniscus-curvature feedback closure (documented approximation -- see
# module docstring). Locally higher pressure -> smaller meniscus radius of
# curvature R = sigma/P -> larger protrusion depth delta -> reduced slip.
# MENISCUS_COEFF is deliberately kept modest (a fraction of 1) so the
# correction acts as a partial, physically-motivated damping of slip rather
# than driving it all the way to zero -- true zero-slip collapse would mean
# the meniscus has depinned/flooded the groove entirely, a different regime
# this simple closure is not meant to describe.
# ---------------------------------------------------------------------------
MENISCUS_COEFF = 0.2   # O(1) sensitivity coefficient, exposed for calibration
MENISCUS_CAP = 0.7     # maximum fractional reduction of lambda0 (never to zero)


def meniscus_protrusion(P, groove: Groove):
    """Depth (m) the gas-liquid meniscus bows into a groove under local
    liquid pressure P (Pa), via Young-Laplace R = sigma/P and a circular-arc
    sagitta over the shear-free (gas) fraction of the period.
    """
    P = np.maximum(P, 1e-6)  # avoid singular R at P -> 0 (flat meniscus)
    R = SIGMA / P
    return groove.phi ** 2 * groove.pitch / (8.0 * R)


def local_slip_lengths(P, groove: Groove, meniscus: bool):
    """Return (lambda(x), lambda_t(x)) [m], the local velocity- and thermal-
    slip lengths given local film pressure P(x).

    Thermal slip is assumed to equal the velocity slip length at the
    flat-meniscus level and to respond to curvature the same way; this
    mirrors the structural symmetry the paper draws between lambda and
    lambda_t but is, again, a simplifying stand-in for the exact (and
    unavailable) thermal-slip closure.
    """
    lam0 = groove.lambda0
    if not meniscus or lam0 == 0.0:
        return np.full_like(P, lam0), np.full_like(P, lam0)
    delta = meniscus_protrusion(P, groove)
    reduction = np.clip(MENISCUS_COEFF * delta / groove.pitch, 0.0, MENISCUS_CAP)
    lam = lam0 * (1.0 - reduction)
    return lam, lam.copy()


# ---------------------------------------------------------------------------
# Slip-modified lubrication BVP for the film pressure P(x)
# ---------------------------------------------------------------------------
def solve_pressure_field(h, H_current, L, dT, groove: Groove, meniscus: bool,
                          n=81, picard_tol=1e-6, picard_maxit=30):
    """Solve d/dx[ K(x) dP/dx ] = -S(x) on x in [-L/2, L/2], P(+-L/2) = 0,
    for a given (spatially uniform) film thickness h [m], with fixed-point
    (Picard) iteration on the meniscus-curvature slip correction.

    K(x) = (h^3 / 12 mu_l) * (1 + 6 lambda(x)/h)   -- slip-modified Poiseuille
                                                       flow conductance
    S(x) = k_l * dT / [ (h + lambda_t(x)) * LATENT_HEAT * RHO_L ]
                                                    -- local melt volume source

    Returns x, P(x), lambda(x), lambda_t(x).
    """
    x = np.linspace(-L / 2.0, L / 2.0, n)
    dx = x[1] - x[0]
    lam = np.full(n, groove.lambda0)
    lamt = np.full(n, groove.lambda0)

    n_i = n - 2  # interior (unknown) nodes; P = 0 imposed at both ends

    for _ in range(picard_maxit):
        K = (h ** 3 / (12.0 * MU_L)) * (1.0 + 6.0 * lam / h)
        S = K_L * dT / ((h + lamt) * LATENT_HEAT * RHO_L)

        K_half = 0.5 * (K[:-1] + K[1:])  # K at i+1/2 midpoints, len n-1 = n_i+1

        # Conservative finite-volume tridiagonal system for interior node i
        # (global index i = m+1, m = 0..n_i-1):
        #   K_half[m]*(P_{i-1}) - (K_half[m]+K_half[m+1])*P_i + K_half[m+1]*P_{i+1}
        #       = -S_i * dx^2
        # with P_0 = P_{n-1} = 0 (Dirichlet), so the boundary sub-/super-diagonal
        # entries simply drop out of the system.
        sub_full = K_half[:n_i] / dx ** 2      # len n_i, sub_full[0] is a BC term (dropped)
        sup_full = K_half[1:n_i + 1] / dx ** 2  # len n_i, sup_full[-1] is a BC term (dropped)
        diag = -(sub_full + sup_full)
        rhs = -S[1:-1]

        ab = np.zeros((3, n_i))
        ab[0, 1:] = sup_full[:-1]   # super-diagonal
        ab[1, :] = diag             # main diagonal
        ab[2, :-1] = sub_full[1:]   # sub-diagonal
        P_interior = solve_banded((1, 1), ab, rhs)

        P = np.zeros(n)
        P[1:-1] = np.maximum(P_interior, 0.0)  # physical: film over-pressure

        lam_new, lamt_new = local_slip_lengths(P, groove, meniscus)
        converged = np.max(np.abs(lam_new - lam)) < picard_tol * max(groove.lambda0, 1e-12)
        lam, lamt = lam_new, lamt_new  # undamped update: converges in ~10-15 iters here
        if converged:
            break

    return x, P, lam, lamt


def force_residual(h, H_current, L, dT, groove: Groove, meniscus: bool):
    x, P, lam, lamt = solve_pressure_field(h, H_current, L, dT, groove, meniscus)
    support = np.trapezoid(P, x)                  # N/m (per unit width)
    target = RHO_S * GRAVITY * H_current * L  # N/m required to hold up solid
    return support - target


def solve_film_thickness(H_current, L, dT, groove: Groove, meniscus: bool,
                          h_lo=1e-7, h_hi=2e-3, h_guess=None):
    """Root-find the (spatially averaged) film thickness h such that the
    integrated film pressure exactly supports the current solid weight.
    Larger h -> lower pressure everywhere, so the residual is monotonically
    decreasing in h and a simple bracket + Brent solve is robust.

    If ``h_guess`` is given (e.g. the film thickness from the previous,
    nearby time step), a tight bracket around it is tried first -- this is
    purely a performance optimization (fewer force_residual evaluations)
    and does not change the solution.
    """
    resid = lambda hh: force_residual(hh, H_current, L, dT, groove, meniscus)

    if h_guess is not None:
        lo, hi = h_guess / 4.0, h_guess * 4.0
        f_lo, f_hi = resid(lo), resid(hi)
        if f_lo >= 0 and f_hi <= 0:
            return brentq(resid, lo, hi, xtol=1e-10, rtol=1e-8)

    f_lo = force_residual(h_lo, H_current, L, dT, groove, meniscus)
    f_hi = force_residual(h_hi, H_current, L, dT, groove, meniscus)
    # Expand bracket if needed (thin/thick film edge cases).
    while f_lo < 0 and h_lo > 1e-9:
        h_lo /= 3.0
        f_lo = force_residual(h_lo, H_current, L, dT, groove, meniscus)
    while f_hi > 0 and h_hi < 1e-1:
        h_hi *= 3.0
        f_hi = force_residual(h_hi, H_current, L, dT, groove, meniscus)
    h = brentq(resid, h_lo, h_hi, xtol=1e-10, rtol=1e-8)
    return h


def melting_rate(h, L, dT, groove: Groove, H_current, meniscus: bool):
    """Average dH*/dt* [m/s] of the remaining solid height, from the
    Stefan condition integrated across the block length."""
    x, P, lam, lamt = solve_pressure_field(h, H_current, L, dT, groove, meniscus)
    q_flux = K_L * dT / (h + lamt)             # local heat flux, W/m^2
    q_avg = np.trapezoid(q_flux, x) / L
    return q_avg / (RHO_S * LATENT_HEAT), (x, P, lam, lamt)


@dataclass
class SimResult:
    t: np.ndarray       # s
    H: np.ndarray       # remaining solid height, m
    h_film: np.ndarray  # film thickness, m
    q_avg: np.ndarray   # average heat flux, W/m^2
    groove: Groove
    meniscus: bool
    dT: float
    L: float
    H0: float


def run_ccm(groove: Groove, dT, L=0.02, H0=0.06, meniscus=True,
            H_stop_frac=0.03, max_steps=4000):
    """Time-march the lumped CCM model until the solid height falls to
    H_stop_frac * H0 (numerical stand-in for full melt-out).
    """
    H = H0
    t = 0.0
    ts, Hs, hs, qs = [t], [H], [], []

    h = solve_film_thickness(H, L, dT, groove, meniscus)
    rate, _ = melting_rate(h, L, dT, groove, H, meniscus)
    hs.append(h)
    qs.append(rate * RHO_S * LATENT_HEAT)

    steps = 0
    while H > H_stop_frac * H0 and steps < max_steps:
        # Adaptive step: no more than ~1% of remaining height per step.
        dt = float(np.clip(0.01 * H / max(rate, 1e-30), 1e-3, 8.0))
        H_new = max(H - rate * dt, 0.0)
        H_eval = max(H_new, 1e-6 * H0)
        h = solve_film_thickness(H_eval, L, dT, groove, meniscus, h_guess=h)
        rate_new, _ = melting_rate(h, L, dT, groove, H_eval, meniscus)
        q_avg = rate_new * RHO_S * LATENT_HEAT

        t += dt
        H = H_new
        rate = rate_new
        ts.append(t)
        Hs.append(H)
        hs.append(h)
        qs.append(q_avg)
        steps += 1

    return SimResult(t=np.array(ts), H=np.array(Hs), h_film=np.array(hs),
                      q_avg=np.array(qs), groove=groove, meniscus=meniscus,
                      dT=dT, L=L, H0=H0)


def stefan_number(dT):
    return CP_L * dT / LATENT_HEAT


def dT_from_stefan(St):
    return St * LATENT_HEAT / CP_L
