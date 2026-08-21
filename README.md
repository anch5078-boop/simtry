# simtry — CCM Simulation

A numerical simulation of **close-contact melting (CCM) on structured
slip surfaces**, inspired by:

> S. Hu, N. Hu, Y. Lai, Z. Li, X. Gao, L. Fan, "Close-contact melting
> regulated by structured slip surfaces," *Appl. Phys. Lett.* **128**,
> 073903 (2026). https://doi.org/10.1063/5.0311034

## What the paper does

A solid PCM block (ice) melts under its own weight on a heated, grooved
superhydrophobic plate. A thin liquid film forms between the unmelted
solid and the plate; its thickness self-adjusts so the film's pressure
supports the solid's weight (a Stefan/lubrication problem). The grooves
give the film an effective **slip length**, but the gas-liquid meniscus
trapped in each groove bows under the film pressure (Young–Laplace),
which changes the local slip — a "refined" model with this
meniscus-curvature feedback predicts systematically *lower* melting
rates than a "conventional" model that treats the slip length as a
constant, geometry-only property.

## What's implemented here, and why it isn't a byte-for-byte reproduction

The paper's own reduced-order model (its Eqs. 2–3) is a nondimensional
system whose closure — the exact formulas for the flat-meniscus slip
length and its curvature correction, for a given groove period and gas
fraction — lives in the paper's Supplementary Material (Secs. S1–S8) and
an earlier companion paper, **neither of which was available** (only the
7-page main article was supplied to build this).

Rather than guess at unpublished coefficients, `ccm_slip/model.py`
re-derives an equivalent model from first principles:

1. **Stefan (energy) balance** at the melting front, with an effective
   thermal slip length `λt` acting as an added thermal resistance
   (`q ∝ ΔT/(h+λt)`) — bigger `λt` slows melting, matching the paper's
   structure `dH/dτ ∝ -1/(Λ+λt)`.
2. **Slip-modified 1-D lubrication (Reynolds) equation** for the film
   pressure `P(x,t)` along the drainage direction — this produces the
   *spatial* pressure variation (and hence spatially varying meniscus
   curvature) that is the paper's central point, rather than a single
   lumped pressure value.
3. **Force balance**: `∫P dx` supports the current weight of the
   unmelted solid, self-consistently setting the film thickness at each
   instant (root-solved every step).
4. **Meniscus-curvature feedback**: locally higher pressure bows the
   meniscus further into the groove (`R = σ/P`, Young–Laplace), which
   *reduces* the local slip length relative to its flat-meniscus value.

The flat-meniscus slip length itself uses the well-established,
independently citable **Philip (1972) / Lauga & Stone (2003)** closed
form for longitudinal micro-ridges:

```
λ0 = (pitch/π) · ln( sec(π·φ/2) )
```

The curvature-feedback closure (how strongly pressure reduces slip) is
the one genuinely approximate, tunable piece (`MENISCUS_COEFF`,
`MENISCUS_CAP` in `model.py`) — it stands in for the paper's
unavailable perturbation-theory result, capped so it damps slip rather
than driving it to a nonphysical total collapse.

**Bottom line:** this reproduces the paper's *mechanism* and its
qualitative claims (grooved slip surfaces can enhance, leave unchanged,
or suppress CCM melting depending on geometry; the meniscus-resolving
model melts slower than the constant-slip model) using a self-consistent,
independently-derived physical model — not the paper's exact numbers.

## Model self-check: why some surfaces enhance and others suppress

Solving the lumped force balance analytically (uniform film, constant
slip) shows the melting-rate ratio relative to a smooth surface depends
only on `x = h/λ`, the film thickness relative to the slip length:

- `x → ∞` (λ → 0): ratio → 1 (no effect, as expected)
- `x ≈ 2–10` (slip a modest fraction of the film thickness): ratio
  peaks modestly **above** 1 — mild enhancement (few %)
- `x → 0` (λ ≫ h): ratio → 0 — strong **suppression** (thermal slip
  resistance dominates)

`run_paper_like.py` picks three showcase groove geometries via
`pitch_for_lambda0()` that land in each regime, rather than trying to
match the paper's specific Surface A/B/C (whose closure differs from
ours and isn't available). It also demonstrates a genuine model
prediction worth noting: the meniscus feedback can move a surface
*across* this ratio peak — e.g. a surface that's suppressed by the
conventional model can end up *enhanced* by the refined model, because
reducing the (already too-large) slip length pulls it back toward the
peak. The paper's own text only reports its refined model giving
*generally* lower rates than its conventional model; it doesn't rule
out this kind of regime crossing for other geometries.

## Running it

```bash
pip install -r requirements.txt
python3 -m ccm_slip.run_paper_like
```

Produces (in `ccm_slip/figures/`):

- `melting_fraction.png` — `H/H0` vs. normalized time, smooth vs. three
  showcase grooved surfaces, conventional vs. refined model (analogous
  to the paper's Fig. 3a–b)
- `heat_flux.png` — average film heat flux vs. time (analogous to
  Fig. 3c–d)
- `summary.txt` — melting-rate % change relative to the smooth surface,
  for both models, at both Stefan numbers

A single run (`ccm_slip/model.py`) takes ~5–20 s; the full batch script
(16 runs) takes ~35–50 s.

## Files

- `ccm_slip/model.py` — physics: `Groove` geometry, slip-length
  closures, the pressure BVP solver (finite-volume + Picard fixed-point
  iteration on the meniscus feedback), the force-balance root-find, and
  the time-marching `run_ccm()`.
- `ccm_slip/run_paper_like.py` — driver that runs the parameter sweep
  and produces the figures/summary described above.

## Extending this with the real paper closure

If you obtain the paper's Supplementary Material (Secs. S1–S8) or its
Ref. 23 (Hu, Fan, Gao, Hu, *J. Fluid Mech.* **1010**, A46, 2025), the
exact `λ⁽⁰⁾`, `λ⁽¹⁾`, `λt⁽⁰⁾`, `λt⁽¹⁾` formulas can be dropped in to
replace `Groove.lambda0` and `local_slip_lengths()` in `model.py`
without touching the rest of the solver (energy balance, pressure BVP,
force balance, time-marching are all closure-agnostic).
