"""
validate_conservation.py -- energy-conservation check for Model B
(radial_enthalpy_model.py), in the same spirit as
PCM_MATLAB/validateStefan.m's validation of the Version-1 solver: this
codebase should hold itself to the same standard.

What it checks: with both boundaries adiabatic during a pulse OFF
interval (no source, no sink anywhere in the domain), the total
sensible+latent enthalpy

    H(t) = sum_i [ rho_s*cp*T_i(t) + rho_s*Lf*fl_i(T_i(t)) ] * V_i

must stay constant between two time points that are both within the
same OFF interval, to within time-discretization/Picard-convergence
error. This directly exercises the exact defect this script was
written to catch: an earlier version of simulate() had a sign error in
how the conductance matrix's off-diagonal terms were passed to its
tridiagonal solver, which still produced a *converged-looking* Picard
iteration but silently violated conservation -- H measurably (not just
to floating point) drained away during every OFF interval, with no
error or warning anywhere else. That failure mode would not have been
caught by eyeballing T(r,t) plots alone.

Usage:
    python3 validate_conservation.py
Exits with a nonzero status (and prints FAIL) if the check does not
pass, like validateStefan.m's pass/fail report.
"""

import dataclasses
import sys
import numpy as np

from parameters import Params
import radial_enthalpy_model as modelB
from radial_enthalpy_model import liquid_fraction

REL_TOL = 0.02  # max allowed |dH| during an OFF step, relative to the
                # typical |dH| magnitude seen during an ON step in the
                # same run (a "how much does this look like a real
                # heat-input step" scale, not an absolute number)


def main():
    # Short run (a few pulse cycles) at fine save resolution -- this
    # check only needs the early transient, not a full charge.
    P = dataclasses.replace(Params(), t_max=200.0)
    b = 45.0e-6
    res = modelB.simulate(P, b, save_every=1)

    dr = res.r_c[1] - res.r_c[0]
    r_face = P.r_h + np.arange(len(res.r_c) + 1) * dr
    V = np.pi * (r_face[1:] ** 2 - r_face[:-1] ** 2)

    fl_hist = np.array([liquid_fraction(T, P.Ts, P.Tl)[0] for T in res.T_history])
    H = np.sum((P.rho_s * P.cp * res.T_history + P.rho_s * P.Lf * fl_hist) * V[None, :], axis=1)
    # save_every=1 here, so t_saved/T_history line up 1:1 with the per-step
    # pulse_on record. dH[i] = H[i+1]-H[i] is produced by the solve that
    # generated T_history[i+1], whose boundary condition is pulse_on[i+1]
    # -- so index with pulse_on[1:], not pulse_on[:-1]/pulse_state(t_saved),
    # either of which would mislabel every step spanning an ON/OFF transition.
    is_on = res.pulse_on[1:]
    dH = np.diff(H)

    if not np.any(is_on) or not np.any(~is_on):
        print("FAIL: run did not include both an ON and an OFF step -- "
              "check parameters.t_on/t_off/t_max.")
        sys.exit(1)

    scale = np.mean(np.abs(dH[is_on]))
    off_drift = np.max(np.abs(dH[~is_on]))
    rel_drift = off_drift / scale

    print(f"Typical |dH| per step while ON  : {scale:10.4f}  (energy actually entering)")
    print(f"Worst   |dH| per step while OFF : {off_drift:10.4f}  (should be ~0 -- adiabatic)")
    print(f"Ratio (worst OFF drift / typical ON step): {rel_drift:.4%}")

    if rel_drift < REL_TOL:
        print(f"PASS: OFF-interval energy drift is {rel_drift:.4%} of a typical ON-step "
              f"input, well under the {REL_TOL:.0%} tolerance -- Model B conserves energy.")
    else:
        print(f"FAIL: OFF-interval energy drift is {rel_drift:.4%} of a typical ON-step "
              f"input, exceeding the {REL_TOL:.0%} tolerance.")
        sys.exit(1)


if __name__ == "__main__":
    main()
