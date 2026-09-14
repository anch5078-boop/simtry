# PCM_MATLAB_radial — hot pipe / PCM / pulsed copper sleeve / PCM / insulated wall

MATLAB/Octave port of `PCM_python/` (this repo's radial 1-D conduction +
phase-change model for the task diagram: a 1" hot-water pipe at the
centre, PCM filling the annulus out to an insulated container wall,
with a thin copper sleeve embedded partway across the gap and
pulse-heated to help melt the PCM faster). It reuses the exact
numerical scheme already validated in `PCM_MATLAB/` (smoothed liquid
fraction, **energy-conserving secant apparent heat capacity**,
harmonic-mean face conductivity, implicit backward-Euler + Picard
iteration — see `PCM_MATLAB/README.md` for why the secant form matters)
on a 1-D cylindrical grid, since this geometry is axisymmetric and
uniform along the pipe axis.

**Cross-checked against `PCM_python/`**: run with the same parameters,
this MATLAB/Octave port reproduces the Python model's numbers exactly
(e.g. the as-drawn case's `Pmin = 290.6 W/m`, `t50 = 3586 s`,
`t99 = 10254 s`, `t999 = 10426 s` match to the last digit) — the two
are independent implementations of the same discretization, not a
copy-paste, so agreement is a genuine cross-validation of the model.

Runs unmodified in MATLAB or GNU Octave (tested on Octave 8.4, same
version as `PCM_MATLAB/`).

## Quick start

**One file, everything, no dependencies** — `PCM_master_simulation.m` is
the whole model (every function plus both drivers) in a single
self-contained file. Copy just this one file anywhere and run it:

```
octave-cli PCM_master_simulation.m
```

(or open it in MATLAB and run it directly). It runs Part 1 (the
geometry as drawn) then Part 2 (the <3-minute redesign) and writes
everything to `./results_master/asDrawn/` and `./results_master/fastMelt/`.
Takes a few minutes end to end (Part 2's gap-width × power sweep is
the slow part). To change a case, edit `parametersRadial()` near the
top of the file (Part 1) or the `TARGET_S` / `T_SAFETY_CAP` /
`REALISTIC_POWER_CAP` constants where Part 2's driver code begins —
nothing else needs to change.

**Or, the original multi-file layout** (same model, split into one
function per file, `results/` and `results_fast_melt/` instead of
`results_master/...`) if you'd rather read or edit it a function at a
time:

```
cd PCM_MATLAB_radial
octave-cli mainAsDrawn.m     % as-drawn case + reference + position sweep -> results/
octave-cli mainFastMelt.m    % re-dimensioned <3-minute design -> results_fast_melt/
```

(or open either script in MATLAB and run it directly). Both layouts
are kept in sync and produce identical numbers — `PCM_master_simulation.m`
was assembled directly from these files.

## The two cases

### `mainAsDrawn.m` — the geometry as drawn

Sleeve centered in a 32.3 mm pipe→wall gap, pulse-heated at the lowest
power that reaches the 40 °C melt point within 300 s of a cold start,
then bang-bang-held there. Also runs a no-sleeve reference and a
sleeve-position sweep to find the minimum achievable melting time for
this architecture.

| Case | sleeve position | min. sleeve peak power | time to fully melt |
|---|---|---|---|
| No sleeve (pipe conduction only) | — | — | 3 h 33 m 36 s |
| As drawn (centered) | 50% of gap | 291 W/m | 2 h 53 m 46 s |
| **Minimum melting time found** | **85% of gap** | 375 W/m | **2 h 20 m 04 s** |

See `PCM_python/README.md` (same numbers, same discussion) for why the
optimum sits near the outer wall rather than in the middle: the PCM's
low conductivity means conduction time scales with distance², so the
sleeve helps most where it shortens the *longest* remaining diffusion
path — the outer region, farthest from the hot pipe.

Output (`results/`): `fig1_liquid_fraction.png`, `fig2_position_sweep.png`,
`fig3_temperature_profiles.png`, `fig4_sleeve_pulsing.png`, `sweep.csv`,
`summary.txt`.

### `mainFastMelt.m` — re-dimensioned for melting in under 3 minutes

A follow-up design question: keeping the pipe radius and PCM material
fixed, how small does the annulus need to be, and what sleeve power
does it take, to melt the *whole* charge in under 180 s? Conduction
time through this low-k PCM scales with distance², so this is a real
re-dimensioning (millimetres, not tens of millimetres), not a tweak.

The script:
1. sweeps gap width with **no sleeve** (pure pipe conduction) — needs
   the gap down to ~4 mm to melt unaided in <180 s;
2. sweeps gap width **with** the sleeve (re-optimizing its position for
   this faster regime — the optimum moves to ~70% of the gap, not the
   ~85% found above), bisecting the minimum sleeve peak power that
   still clears the target at each gap;
3. picks the largest gap whose required power stays under a realistic
   2 kW/m cap, adds 5% margin, and verifies on a finer grid.

**Recommended design:**

| Quantity | Value |
|---|---|
| Pipe radius (fixed, given) | 12.7 mm |
| **PCM annulus gap width** | **12 mm** (was 32.3 mm) |
| Container inner radius | 24.7 mm (was 45 mm) |
| Copper sleeve | 1.0 mm thick, at 70% of the gap |
| Sleeve pulse power | **1700 W/m** (≈510 W over 0.3 m), 180 °C safety cutoff |
| **Result** | **fully melted in ~167 s (2.8 min)** |

At this same 12 mm gap, pipe conduction *alone* reaches only ~40%
melted after 400 s — the sleeve is what makes the 3-minute target
reachable at all while keeping a non-trivial PCM layer. Unlike the
as-drawn case, the sleeve here runs at essentially **continuous full
power** for the whole transient rather than duty-cycling on/off (it
never reaches its 180 °C safety cutoff within 167 s) — the objective
changed from standby efficiency to raw speed. Cheaper/slower
alternatives (e.g. ~10 mm gap at ~1.0 kW/m, ~8 mm at ~0.6 kW/m) are in
`results_fast_melt/sleeve_sweep.csv`.

Output (`results_fast_melt/`): `fig1_gap_sweeps.png` (the two tradeoff
curves), `fig2_liquid_fraction.png`, `fig3_temperature_profiles.png`
(two fronts meeting — from the pipe and from the sleeve),
`fig4_sleeve_pulsing.png`, `no_sleeve_sweep.csv`, `sleeve_sweep.csv`,
`summary.txt`.

## File structure

| File | Role |
|---|---|
| `PCM_master_simulation.m` | **Everything in one file** — every function below plus both drivers, self-contained (see "Quick start") |
| `parametersRadial.m` | All physical/geometric/control constants (edit this to change the as-drawn case) |
| `materialStruct.m` | Packs P into the small struct solveOneStepRadial.m needs (mushy bounds, effective rho·cp / rho·L) |
| `sleeveRadii.m` | Inner/outer sleeve radius for a given fractional position across the gap |
| `createGeometryRadial.m` / `createGeometryNoSleeve.m` | Uniform radial finite-volume grid + sleeve cell mask (or no sleeve at all) |
| `liquidFractionRadial.m`, `liquidFractionDerivRadial.m` | Smoothed liquid fraction and its derivative |
| `apparentCapacityRadial.m` | Energy-conserving secant apparent heat capacity |
| `solveOneStepRadial.m` | One implicit conduction/phase-change time step (Dirichlet pipe BC, adiabatic wall BC, sleeve heat source) |
| `volumeAvg.m` | Volume-weighted average over a cell mask (cell volumes vary with radius on this grid) |
| `rampFinalSleeveTemp.m`, `findMinPowerRadial.m` | Bisect for the minimum sleeve peak power that reaches Tm within a ramp-time target |
| `runFullMeltRadial.m` | Full melt time-marching loop with the bang-bang thermostat (efficiency *or* speed framing via `Tsetpoint`) |
| `runPositionSweepRadial.m` | Sleeve-position sweep for the as-drawn case |
| `t999ForPower.m`, `minPowerForTarget.m`, `positionSweepQuick.m` | Fast-melt-specific helpers: bisect for min. power to hit a *time* target, and rank sleeve positions |
| `sweepNoSleeveGaps.m`, `sweepWithSleeveGaps.m` | Gap-width sweeps for the fast-melt design |
| `mainAsDrawn.m` + `makePlotsAsDrawn.m` | Driver + plots for the as-drawn case |
| `mainFastMelt.m` + `makePlotsFastMelt.m` | Driver + plots for the <3-minute redesign |
| `numOrEmpty.m`, `formatTimeOrNARadial.m` | Small formatting utilities (kept as separate files, not script-local functions, for MATLAB/Octave compatibility — same reason as `PCM_MATLAB/formatTimeOrNA.m` et al.) |

## Assumptions

Same as `PCM_python/params.py` (container radius, sleeve thickness,
PCM dataset, initial temperature, ramp-up target) — see that file's
table for the full list and reasoning. Edit `parametersRadial.m` (as-drawn
case) or the constants at the top of `mainFastMelt.m` (fast-melt case)
to try different numbers; nothing else needs to change.
