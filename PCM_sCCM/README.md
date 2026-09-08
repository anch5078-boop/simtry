# PCM_sCCM — Version 2: slip-enhanced close-contact melting (vertical tube, pulsed heater)

Python companion to `PCM_MATLAB` (Version 1, conduction-only) modeling a
**vertical-tube** PCM thermal battery charged by a **pulsed** cylindrical
heater sleeve, where melting is dominated by **slip-enhanced close-contact
melting (sCCM)** rather than plain conduction.

Grounding reference: Li et al. 2026, *"Pulse heating and slip enhance
charging of phase-change thermal batteries"* (Nature) — the companion
paper this repo's CCM work is being extended to track. That paper
reports a slip coating giving a Navier slip length **b ≈ 45–90 µm** via
sCCM, used here directly as the slip-length range (`Params.slip_scenarios`
in `parameters.py`).

## Why this needed a new model, not just an extension of PCM_MATLAB

`PCM_MATLAB` (Version 1) is a 3D **conduction-only** solver for a PCM
sitting around two heater sleeves — no gravity-driven melt motion at
all. That is the right model for a flat/gravity-drainage close-contact
melting setup, but not for this one: here the tube is **vertical**, so
gravity acts *along* the axis the heater runs up, not across the melt
film. That changes the governing mechanism entirely — see
`ccm_slip_model.py`'s module docstring for the physical picture (a solid
PCM column sinking under gravity, lubricated by a thin melt film against
the heater surface) and its full derivation.

## The two models

| | Model A — `ccm_slip_model.py` | Model B — `radial_enthalpy_model.py` |
|---|---|---|
| Type | 0-D lumped ODE (algebraic force balance + mass/energy balance) | 1-D transient radial enthalpy PDE (finite volume, implicit) |
| Geometry | Solid column **sinks** past the heater as it melts (matches the real device mechanism) | Fixed annulus, melt front grows radially outward, **no sinking** |
| Gives you | Overall charge-time trends: melt fraction, t50/t90/t95/t99, film thickness, heat flux, sinking velocity | Actual T(r,t) fields and melt-front r(t), including pulsing detail |
| Slip enters via | A rigorously re-derived one-wall-Navier-slip lubrication force balance (closed form, reduces to the classic no-slip CCM 1/4-power scaling at b=0 — see docstring) | An engineering effective-conductivity blend (slip-decay term + Raithby-Hollands-type natural-convection correlation) — explicitly flagged as an approximation, not a first-principles derivation |
| Run time | ~2 s for all 3 slip scenarios | ~5 s for all 3 slip scenarios |

Run both and compare:

```bash
cd PCM_sCCM
pip install -r requirements.txt   # numpy, scipy, matplotlib
python3 run_all.py
```

This prints a t50/t90/t95/t99 milestone table for each model and slip
scenario, and writes six comparison plots to `results/`.

To re-check Model B's energy conservation (see below) after changing
the mesh, time step, or Picard settings:

```bash
python3 validate_conservation.py
```

## Key result from the shipped placeholder parameters

With the placeholder geometry/properties in `parameters.py` (see below
before trusting the numbers quantitatively; exact numbers reproduce via
`python3 run_all.py`):

|  | t50 | t90 | t95 |
|---|---|---|---|
| Model A, no slip | 52.4 min | 109.5 min | not reached |
| Model A, b=45 µm | 48.5 min | 102.2 min | 114.1 min |
| Model A, b=90 µm | 46.3 min | 97.2 min | 108.3 min |
| Model B, no slip | 46.4 min | 66.3 min | 69.0 min |
| Model B, b=45 µm | 46.3 min | 66.3 min | 68.8 min |
| Model B, b=90 µm | 46.3 min | 66.3 min | 68.5 min |

**Do not read the Model A vs. Model B *absolute* times as "which
mechanism is faster"** — they are not answering the same question. Model
A charges an entire `H_col` = 150 mm solid column past a 40 mm heated
section; Model B fills only its own fixed `R - r_h` = 24 mm annulus
(per unit height, no axial extent at all). Their total thermal mass and
geometric scope were not chosen to match, so their timescales aren't
comparable — only trends *within* each model are.

What *is* comparable, and is the actual point of running both models:

* **Model A** (the physically-correct mechanism for this geometry): the
  sinking column keeps the film thin (~250–350 µm) and the heat flux
  high (~14–20 kW/m²) for the *entire* charge, because sinking
  continuously renews the thin film instead of letting it thicken. Slip
  gives a modest but **sustained, monotonically-growing-with-b** benefit
  all the way to complete melt-out: a **t90 speed-up factor of
  1.07–1.13x** (i.e. ~7–13% faster) across the paper's b=45–90 µm range
  (`A_melt_fraction_vs_time.png`, `A_film_and_flux.png`; exact ratios
  printed by `run_all.py`).
* **Model B** (fixed-position radial melting — no sinking): slip's
  effect is real but **transient** — clearly visible for the first ~5
  minutes while the gap is still comparable to b
  (`B_melt_front_pulsing.png`, where b=90 µm is consistently ahead of
  b=0) but essentially washed out by the time t90 is reached, once the
  gap has grown well past the slip-length scale and natural convection
  dominates instead (`B_liquid_fraction_vs_time.png` — the three curves
  overlap almost exactly).
* This contrast **is** the point: it is a quantitative illustration of
  why the sinking-solid mechanism (Model A / the reference paper's
  actual device) keeps extracting a benefit from the same slip coating
  for the whole charge, where a static radial-melting geometry would
  only benefit from it briefly, at the very start.
* `B_melt_front_pulsing.png` also shows a real, physically-legitimate
  (energy-conservation-checked — see below) effect worth knowing about
  when designing a pulse schedule: during each OFF interval the melt
  front does not just stall, it **partially retreats** — stored
  superheat in the near-wall liquid diffuses outward and forward
  instead of being preserved, so the front growth curve is a sawtooth,
  not a staircase. Very long OFF periods relative to the melt-layer's
  thermal diffusion time trade away some of what the ON period gained.

## PLACEHOLDER DATA — replace before trusting numbers quantitatively

Every geometry and material constant lives in one place, `parameters.py`
(`Params` dataclass), mirroring `PCM_MATLAB/parameters.m`'s
single-source-of-truth convention. They are realistic-magnitude
placeholders (same paraffin-type PCM / copper heater sleeve as the
Version-1 default), **not measured values**. Two properties new to this
model (not needed by the conduction-only Version 1) matter most and are
flagged in `parameters.py`:

* `mu_l` (liquid viscosity) — sets the CCM film thickness (`delta ~
  mu_l**0.25`, so this is the single most consequential number you
  should replace with a real datasheet/measured value).
* `beta_l` (liquid thermal expansion coefficient) — sets the
  natural-convection Rayleigh number in Model B.

Swap in your real tube radius, container radius, column height,
heater temperature/pulse schedule, and PCM thermophysical properties in
`parameters.py` — nothing else in the code should need to change.

## Modeling assumptions (stated explicitly — read before citing numbers)

Both modules' docstrings spell these out in full; the headline ones:

* **Model A** lumps dT and film thickness as spatially uniform over the
  wetted heater height at each instant (only the axial *flow rate*
  varies) — a standard first simplification in the CCM literature.
  Pulsing is modeled by neglecting the film's own thermal mass (OFF ⇒
  q=0 instantly, no "coasting"). The force balance is undefined once
  the driving weight vanishes near complete melt-out, so t99 is reported
  but flagged as unreliable for this model (t95 is the more robust
  summary number) — see `ccm_slip_model.py`.
* **Model B**'s slip-plus-convection effective-conductivity submodel is
  an explicitly-labeled engineering approximation (not a rigorous
  derivation like Model A's force balance) standing in for physics
  (thin-film lubrication flow, buoyant circulation) that a 1-D model
  cannot resolve directly — see `radial_enthalpy_model.py`.
* Both models use one fixed density for capacity bookkeeping (no
  melting volume-change/expansion physics) and an adiabatic outer
  container wall, consistent with `PCM_MATLAB`'s Version-1 conventions.
* Model B's numerics were checked directly against energy conservation
  (total sensible+latent enthalpy tracked through a full adiabatic OFF
  window) during development, the same standard `PCM_MATLAB` holds
  itself to via `validateStefan.m` — worth re-checking if you change the
  Picard iteration count/tolerance or the mesh.

## File structure

| File | Role |
|---|---|
| `parameters.py` | All physical/geometric/numerical constants (`Params` dataclass) — edit this to change a case |
| `ccm_slip_model.py` | Model A: 0-D lumped sinking-column force-balance sCCM model, with full derivation in the module docstring |
| `radial_enthalpy_model.py` | Model B: 1-D radial transient enthalpy-method PDE solver (secant apparent-capacity method, same energy-conserving technique as `PCM_MATLAB/apparentCapacity.m`), pulsed Dirichlet/adiabatic wall BC, slip+convection effective-conductivity submodel |
| `run_all.py` | Driver: runs both models across the three slip scenarios, prints the milestone table, saves `results/*.png` |
| `validate_conservation.py` | Energy-conservation check for Model B (same role as `PCM_MATLAB/validateStefan.m`) |
| `requirements.txt` | Python dependencies (`numpy`, `scipy`, `matplotlib`) |
| `results/` | Generated comparison plots (not hand-edited; regenerate with `run_all.py`) |
