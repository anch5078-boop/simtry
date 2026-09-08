# PCM_sCCM_MATLAB — Version 2: slip-enhanced close-contact melting (vertical tube, pulsed heater)

MATLAB/Octave port of `PCM_sCCM` (the Python version of this Version-2
model), written in the same per-function-per-file style as
`PCM_MATLAB` (Version 1, conduction-only). Models a **vertical-tube**
PCM thermal battery charged by a **pulsed** cylindrical heater sleeve,
where melting is dominated by **slip-enhanced close-contact melting
(sCCM)** rather than plain conduction.

Grounding reference: Li et al. 2026, *"Pulse heating and slip enhance
charging of phase-change thermal batteries"* (Nature) — reports a slip
coating giving a Navier slip length **b ≈ 45–90 µm** via sCCM, used
here directly as the slip-length range (`runAll.m`'s `scenarioB`).

Runs unmodified in MATLAB or GNU Octave (tested on Octave 8.4).

## Quick start

```matlab
cd PCM_sCCM_MATLAB
runAll        % runs both models across 3 slip scenarios, prints the
              % milestone table, saves results/*.png
```

```matlab
validateConservation   % energy-conservation regression check for Model B
```

`runAll` takes **~1 minute** in Octave (Model A's root-find-per-step
loop dominates; MATLAB's JIT will likely be faster). If you only need
one model/scenario:

```matlab
P = parameters();
R = ccmSlipModel(P, 45e-6);        % Model A, b = 45 um
R = radialEnthalpyModel(P, 90e-6); % Model B, b = 90 um
```

## Why this needed a new model, not just an extension of PCM_MATLAB

`PCM_MATLAB` (Version 1) is a 3D **conduction-only** solver for a PCM
sitting around two heater sleeves — no gravity-driven melt motion at
all. That is the right model for a flat/gravity-drainage close-contact
melting setup, but not for this one: here the tube is **vertical**, so
gravity acts *along* the axis the heater runs up, not across the melt
film. That changes the governing mechanism entirely — see
`ccmSlipModel.m`'s header comment for the physical picture (a solid PCM
column sinking under gravity, lubricated by a thin melt film against
the heater surface) and its full derivation.

## The two models

| | Model A — `ccmSlipModel.m` | Model B — `radialEnthalpyModel.m` |
|---|---|---|
| Type | 0-D lumped ODE (algebraic force balance + mass/energy balance) | 1-D transient radial enthalpy PDE (finite volume, implicit) |
| Geometry | Solid column **sinks** past the heater as it melts (matches the real device mechanism) | Fixed annulus, melt front grows radially outward, **no sinking** |
| Gives you | Overall charge-time trends: melt fraction, t50/t90/t95/t99, film thickness, heat flux, sinking velocity | Actual T(r,t) fields and melt-front r(t), including pulsing detail |
| Slip enters via | A rigorously re-derived one-wall-Navier-slip lubrication force balance (closed form, reduces to the classic no-slip CCM 1/4-power scaling at b=0 — see header comment) | An engineering effective-conductivity blend (slip-decay term + Raithby-Hollands-type natural-convection correlation) — explicitly flagged as an approximation, not a first-principles derivation |
| Root-finder | `fzero` (bracketed on `[1e-9, 5e-2]` m) | n/a — direct sparse linear solve each Picard iterate |

## Key result from the shipped placeholder parameters

Reproduces via `runAll`:

|  | t50 | t90 | t95 |
|---|---|---|---|
| Model A, no slip | 52.4 min | 109.5 min | not reached |
| Model A, b=45 µm | 48.5 min | 102.2 min | 114.0 min |
| Model A, b=90 µm | 46.3 min | 97.2 min | 108.3 min |
| Model B, no slip | 46.4 min | 66.3 min | 69.0 min |
| Model B, b=45 µm | 46.3 min | 66.3 min | 68.8 min |
| Model B, b=90 µm | 46.3 min | 66.3 min | 68.5 min |

**Do not read the Model A vs. Model B *absolute* times as "which
mechanism is faster"** — they are not answering the same question.
Model A charges an entire `P.Hcol` = 150 mm solid column past a 40 mm
heated section; Model B fills only its own fixed `R - rH` = 24 mm
annulus (per unit height, no axial extent at all). Their total thermal
mass and geometric scope were not chosen to match, so their timescales
aren't comparable — only trends *within* each model are.

What *is* comparable, and is the actual point of running both models:

* **Model A** (the physically-correct mechanism for this geometry): the
  sinking column keeps the film thin (~250–350 µm) and the heat flux
  high (~14–20 kW/m²) for the *entire* charge, because sinking
  continuously renews the thin film instead of letting it thicken. Slip
  gives a modest but **sustained, monotonically-growing-with-b**
  benefit all the way to complete melt-out: a **t90 speed-up factor of
  1.07–1.13x** across the paper's b=45–90 µm range
  (`results/A_melt_fraction_vs_time.png`, `results/A_film_and_flux.png`;
  exact ratios printed by `runAll`).
* **Model B** (fixed-position radial melting — no sinking): slip's
  effect is real but **transient** — clearly visible for the first ~5
  minutes while the gap is still comparable to b
  (`results/B_melt_front_pulsing.png`, where b=90 µm is consistently
  ahead of b=0) but essentially washed out by the time t90 is reached,
  once the gap has grown well past the slip-length scale and natural
  convection dominates instead (`results/B_liquid_fraction_vs_time.png`
  — the three curves overlap almost exactly).
* This contrast **is** the point: it is a quantitative illustration of
  why the sinking-solid mechanism (Model A / the reference paper's
  actual device) keeps extracting a benefit from the same slip coating
  for the whole charge, where a static radial-melting geometry would
  only benefit from it briefly, at the very start.
* `results/B_melt_front_pulsing.png` also shows a real,
  physically-legitimate (energy-conservation-checked — see below)
  effect worth knowing about when designing a pulse schedule: during
  each OFF interval the melt front does not just stall, it **partially
  retreats** — stored superheat in the near-wall liquid diffuses
  outward and forward instead of being preserved, so the front growth
  curve is a sawtooth, not a staircase. Very long OFF periods relative
  to the melt-layer's thermal diffusion time trade away some of what
  the ON period gained.

## PLACEHOLDER DATA — replace before trusting numbers quantitatively

Every geometry and material constant lives in one place,
`parameters.m` (the `Params` struct), mirroring
`PCM_MATLAB/parameters.m`'s single-source-of-truth convention. They are
realistic-magnitude placeholders (same paraffin-type PCM / copper
heater sleeve as the Version-1 default), **not measured values**. Two
properties new to this model (not needed by the conduction-only
Version 1) matter most and are flagged in `parameters.m`:

* `P.muLiq` (liquid viscosity) — sets the CCM film thickness
  (`delta ~ muLiq^(1/4)`, so this is the single most consequential
  number you should replace with a real datasheet/measured value).
* `P.betaLiq` (liquid thermal expansion coefficient) — sets the
  natural-convection Rayleigh number in Model B.

Swap in your real tube radius, container radius, column height, heater
temperature/pulse schedule, and PCM thermophysical properties in
`parameters.m` — nothing else in the code should need to change.

## Modeling assumptions (stated explicitly — read before citing numbers)

Both `ccmSlipModel.m` and `radialEnthalpyModel.m`'s header comments
spell these out in full; the headline ones:

* **Model A** lumps dT and film thickness as spatially uniform over the
  wetted heater height at each instant (only the axial *flow rate*
  varies) — a standard first simplification in the CCM literature.
  Pulsing is modeled by neglecting the film's own thermal mass (OFF ⇒
  q=0 instantly, no "coasting"). The force balance is undefined once
  the driving weight vanishes near complete melt-out, so t99 is
  reported but flagged as unreliable for this model (t95 is the more
  robust summary number).
* **Model B**'s slip-plus-convection effective-conductivity submodel
  (`enhancementFactor.m`) is an explicitly-labeled engineering
  approximation (not a rigorous derivation like Model A's force
  balance) standing in for physics (thin-film lubrication flow, buoyant
  circulation) that a 1-D model cannot resolve directly.
* Both models use one fixed density for capacity bookkeeping (no
  melting volume-change/expansion physics) and an adiabatic outer
  container wall, consistent with `PCM_MATLAB`'s Version-1 conventions.
* Model B's numerics were checked directly against energy conservation
  (total sensible+latent enthalpy tracked through a full adiabatic OFF
  window) during development, the same standard `PCM_MATLAB` holds
  itself to via `validateStefan.m` — worth re-checking with
  `validateConservation.m` if you change the Picard iteration
  count/tolerance or the mesh.

## Octave/gnuplot rendering notes

*(Same headless-renderer caveat as `PCM_MATLAB/README.md`.)* Tested
with `octave-cli` under the fallback `gnuplot` graphics toolkit (no Qt
available in this environment): plain line/step plots and `fill()`
shading render correctly, but expect a harmless
`using the gnuplot graphics toolkit is discouraged` warning and a
`iconv failed to convert degree sign` warning (axis/title text uses
"deg C" rather than the ° glyph specifically to avoid the latter). Real
MATLAB, or Octave with the `qt` toolkit, needs neither workaround.

The Model A heat-flux and film-thickness series are genuinely
discontinuous (an instantaneous 0 ↔ ON-value switch every pulse edge —
see `ccmSlipModel.m`'s "Pulsing" assumption above); `runAll.m` plots
those with `stairs()`, not `plot()`, so the figure shows the true
square wave instead of a straight-line-interpolated ramp between
samples.

## File structure

| File | Role |
|---|---|
| `parameters.m` | All physical/geometric/numerical constants (`Params` struct) — edit this to change a case |
| `pulseState.m` | `isOn = pulseState(t, P)` — mirrors `PCM_MATLAB/heaterPower.m`'s pulsing convention |
| `liquidFraction.m` | Smoothed liquid fraction `fl(T)` — identical formula to `PCM_MATLAB/liquidFraction.m` |
| `apparentCapacity.m` | Energy-conserving secant apparent heat capacity — identical formula to `PCM_MATLAB/apparentCapacity.m` |
| `frontPosition.m` | Linear-interpolated melt-front crossing of T=Tm — identical convention to `PCM_MATLAB/frontPosition.m` |
| `deltaEffCubed.m` | The slip-modified effective hydrodynamic gap (cubed) used by Model A's force balance |
| `solveDelta.m` | Root-finds Model A's quasi-steady film thickness each step (`fzero`) |
| `ccmSlipModel.m` | **Model A**: 0-D lumped sinking-column force-balance sCCM model, with full derivation in the header comment |
| `enhancementFactor.m` | Model B's slip+convection effective-conductivity submodel |
| `radialEnthalpyModel.m` | **Model B**: 1-D radial transient enthalpy-method PDE solver (sparse-matrix FV assembly, same triplet pattern as `PCM_MATLAB/solveOneTimeStep.m`), pulsed Dirichlet/adiabatic wall BC |
| `milestonesFromSeries.m` | Shared t50/t90/t95/t99 interpolation helper for both models |
| `runAll.m` | Driver: runs both models across the three slip scenarios, prints the milestone table, saves `results/*.png` |
| `printMilestoneTable.m`, `fmtT.m`, `shadeIntervals.m` | Small shared utilities (kept as separate files rather than script-local functions for MATLAB/Octave compatibility, per `PCM_MATLAB/README.md`'s own convention) |
| `validateConservation.m` | Energy-conservation check for Model B (same role as `PCM_MATLAB/validateStefan.m`) |
| `results/` | Generated comparison plots (not hand-edited; regenerate with `runAll`) |
