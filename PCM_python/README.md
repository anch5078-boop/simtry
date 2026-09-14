# PCM_python — hot pipe / PCM / pulsed copper sleeve / PCM / insulated wall

A 1-D radial (cylindrical) transient conduction + phase-change model of
the geometry in the task diagram: a 1" hot-water pipe at the centre,
PCM filling the annulus out to an insulated container wall, with a
thin copper sleeve embedded partway across the gap and pulse-heated to
help melt the PCM faster.

This is a companion model to `PCM_MATLAB/` (the repo's existing
3-D Cartesian two-heater solver) — it reuses the same validated
numerical scheme (smoothed liquid fraction, **energy-conserving secant
apparent heat capacity**, harmonic-mean face conductivity, implicit
backward-Euler + Picard iteration — see `PCM_MATLAB/README.md` for why
the secant form matters) but on a 1-D cylindrical grid, because this
geometry is axisymmetric and uniform along the pipe axis, so the
extra two Cartesian dimensions in `PCM_MATLAB/` buy nothing here and
cost a lot of runtime. Every "case" below (find minimum sleeve power +
run to full melt) takes well under a second; the whole position sweep
takes well under a minute.

## Quick start

```bash
cd PCM_python
python3 main.py       # builds the case(s), runs everything, writes results/
```

Requires `numpy`, `scipy`, `matplotlib` (`pip install numpy scipy matplotlib`).

## The model

**Geometry** (radial cross-section, per unit length of pipe):
`r_pipe` → PCM (inner) → copper sleeve → PCM (outer) → `r_wall`.

**Boundary conditions:**
- Inner boundary (pipe OD, `r_pipe`): **Dirichlet**, fixed at the
  hot-water temperature. This assumes the pipe wall itself has
  negligible thermal resistance compared to the PCM annulus (thin
  metal wall, turbulent well-mixed water → high internal h) — the
  standard simplification for a first screening model.
- Outer boundary (container wall, `r_wall`): **adiabatic**. The
  container is just a vessel, not an active heat sink.

**Copper sleeve control — "pulse heated at lowest power to reach 40 °C":**
Because the domain has a hot Dirichlet inner wall and no heat sink,
*any* nonzero sleeve power (or even zero) eventually raises the sleeve
above the melt point given enough time — conduction from the pipe
alone gets there eventually. So "reach 40 °C" only becomes a
well-posed design constraint once it is tied to a timescale: this
model targets reaching the melt point within `t_ramp_target = 300 s`
of a cold start (an assumed, editable design target — "heat it up
quickly", not "wait for it to happen").

Given that target, `sweep.find_min_power` **bisects for the smallest
sleeve peak power** `P_min` such that a sleeve run continuously at
that power reaches `Tm` within `t_ramp_target`. The full melt
simulation then runs a **bang-bang thermostat**: heater ON at `P_min`
whenever the sleeve's volume-averaged temperature is below `Tm`, OFF
once it reaches `Tm` — i.e. it pulses at exactly the peak power found,
duty-cycled to just barely hold the sleeve at the PCM melting point
without overheating it. That is the literal "pulsed, lowest power,
reaches 40 °C" operating rule from the prompt, and it is what
`fig4_sleeve_pulsing.png` shows: a short burst of full-power pulses to
ramp up, tapering to a low duty cycle once melting starts around the
sleeve (its own latent-heat sink then does most of the work of holding
it near 40 °C).

## Assumptions (not given in the prompt — edit `params.py`)

| Quantity | Value used | Why |
|---|---|---|
| Container inner radius `r_wall` | 45 mm | Not dimensioned in the sketch; picked for a "wider gap" a few pipe-radii wide, per the diagram's proportions. |
| Copper sleeve thickness | 1.5 mm | Thin sleeve, as drawn. |
| Sleeve position (baseline) | centered in the gap | As drawn ("centered in the gap" in the sketch caption). |
| PCM material | paraffin/RT42-class wax: ρ=830 kg/m³, cp=2100 J/(kg·K), k=0.20/0.15 W/(m·K) (solid/liquid), L=170 kJ/kg | Matches "PCM melt temp 40 °C" — RT42-class commercial paraffin melts in this band. **Replace with your actual PCM datasheet values.** |
| Mushy-zone width | 4 °C (38–42 °C) | Numerical regularization only (smooths the latent-heat spike); shrink it and results are unchanged (checked). |
| Initial temperature | 20 °C uniform | Ambient start. |
| Sleeve ramp-up target | 300 s | "Pulse heat it quickly" — see control discussion above. |
| Grid | 240 radial cells | Grid-independence checked: 240 vs 600 cells changes melting-time results by <1%. |

All of the above are fields in `params.py` — nothing else in the code
needs to change to try different numbers.

## Results (with the assumptions above)

| Case | sleeve position | min. sleeve peak power | **time to 99% melted** | **time to fully melt** |
|---|---|---|---|---|
| No sleeve (pipe conduction only) | — | — | 3 h 30 m | **3 h 33 m 36 s** |
| As drawn (sleeve centered in gap) | 50% of gap | 291 W/m | 2 h 51 m | **2 h 53 m 46 s** |
| **Minimum melting time found** | **85% of gap** (i.e. close to the outer wall) | 375 W/m | 2 h 17 m | **2 h 20 m 04 s** |

("W/m" = watts per metre of pipe length; multiply by the actual heated
pipe length for a total wattage, e.g. ×0.5 m → ~188 W for the optimal
case over a 0.5 m section.)

**Headline: the pulsed copper sleeve is worth doing (saves ~40 min /
~19% over pipe conduction alone even placed exactly as drawn), and
relocating it from the middle of the gap to about 85% of the way from
the pipe to the wall (i.e. near-ish the outer wall, not centered) cuts
the melt time further, to about 2 h 20 m — roughly 34% faster than
letting the pipe do it alone, and ~19% faster than the centered
placement in the sketch.**

Why the optimum sits near the wall rather than in the middle: the
PCM's low conductivity (~0.2 W/m·K) means conduction time scales
roughly with the *square* of the distance it has to travel. The pipe
at 120 °C is already a very strong, "free" heat source for the PCM
near it; the outer PCM (near the adiabatic wall) is the slow part of
the problem because it is both farthest from the pipe and has nothing
helping it from the outside. Moving the sleeve away from the pipe (which
doesn't need the help) and toward the wall (which does) directly
shortens the longest, rate-limiting diffusion path. Push it too close
to the wall, though, and there is so little outer PCM left that the
sleeve is "wasted" heating a thin shell while the inner region (now
proportionally huge) again becomes rate-limiting — hence the sweep
curve (`fig2_position_sweep.png`) is a shallow U with a broad, flat
minimum around 80–85% of the gap, not a knife-edge.

Reproducing this or trying other geometries/materials: edit
`params.py` and rerun `python3 main.py`.

## Output files (`results/`)

- `fig1_liquid_fraction.png` — melting progress (avg. liquid fraction
  vs. time), centered vs. optimal vs. no-sleeve.
- `fig2_position_sweep.png` — melting time vs. sleeve radial position;
  shows the minimum.
- `fig3_temperature_profiles.png` — radial temperature snapshots for
  the optimal case, showing the two melt fronts (from the pipe, and
  from the sleeve) growing toward each other.
- `fig4_sleeve_pulsing.png` — sleeve temperature and pulsed power vs.
  time for the as-drawn (centered) case.
- `sweep.csv` — full position-sweep data (every fraction tried, its
  `P_min`, and t50/t90/t99/t99.9).
- `summary.json` — headline numbers in machine-readable form.

## File structure

| File | Role |
|---|---|
| `params.py` | All physical/geometric/control constants (edit this to change a case) |
| `model.py` | Liquid fraction, secant apparent heat capacity, radial FV grid, one implicit time step |
| `simulate.py` | Build a case, bisect for the minimum sleeve peak power, run the full bang-bang-controlled melt simulation |
| `sweep.py` | Sleeve-position sweep; no-sleeve reference case |
| `main.py` | Driver: runs baseline + reference + sweep, writes plots/CSV/summary |
| `fast_melt_design.py` | Re-dimensions the annulus for <3-minute melting (see "Fast-melt design" below) |

## Fast-melt design: melting in under 3 minutes

`fast_melt_design.py` (`python3 fast_melt_design.py`) answers a
follow-up question: **re-dimension the same architecture (pipe + PCM +
pulsed copper sleeve + PCM + insulated wall) so the whole charge melts
in under 180 s**, keeping the pipe radius and PCM material fixed.

**Why this needs a real redesign, not a tweak:** conduction time
through a low-conductivity PCM (k≈0.2 W/m·K) scales with the *square*
of the distance heat has to travel. The original 32.3 mm gap took
2.3–3.5 h to melt; getting under 3 minutes needs a length scale
roughly `sqrt(170/8500) ≈ 1/7` of that, i.e. a few mm to ~1 cm, not
tens of mm. This is a real, physical constraint of any conduction-only
PCM system with this material — not a numerical artifact — so the
"optimization" is mostly about how much annulus thickness (= usable
PCM thermal-storage volume) can be kept while still clearing the
3-minute bar, for a realistic sleeve power.

**Method:**
1. Sweep gap width with **no sleeve at all** (pure pipe conduction) —
   this needs the gap down to ~4 mm to melt in <180 s unaided
   (`fig1_gap_sweeps.png`, left panel).
2. Sweep gap width **with** the pulsed sleeve (position re-optimized
   for this regime — the optimum moves to ~70% of the gap, not the
   ~85% found for the original slow case), bisecting for the minimum
   sleeve peak power that still clears the target at each gap
   (right panel). This quantifies how much *extra* PCM thickness the
   sleeve buys, and at what power cost.
3. Pick the largest gap whose required power stays under a "realistic"
   cap (2 kW/m — a few hundred W over a practical ≤0.5 m pipe section,
   well within a compact resistive cartridge heater or induction
   coil), add 5% power margin, and verify on a finer grid.

**Recommended fast-melt dimensions:**

| Quantity | Value |
|---|---|
| Pipe radius (fixed, given) | 12.7 mm |
| **PCM annulus gap width** | **12 mm** (was 32.3 mm) |
| Container inner radius | 24.7 mm (was 45 mm) |
| Copper sleeve thickness | 1.0 mm, at 70% of the gap (20.6–21.6 mm) |
| Sleeve pulse power | **1700 W/m** of pipe length (≈510 W over 0.3 m) |
| Sleeve safety cutoff | 180 °C (well below paraffin decomposition, ~200–250 °C) |
| **Result** | **fully melted in 167 s (2.8 min)**, vs. the 180 s target |

At this gap, pipe conduction *alone* (no sleeve) only reaches ~40%
melted by 400 s and would take well over 15–20 minutes to finish —
the sleeve is not optional at this size, it is what makes the
3-minute target reachable at all while keeping a non-trivial (12 mm)
PCM layer. Note also the operating character changes from the
"efficiency" case above: to melt this fast, the sleeve runs at
**continuous full power** for essentially the whole transient (it
never reaches its 180 °C safety cutoff within 167 s — see
`fig4_sleeve_pulsing.png`) rather than duty-cycling on/off. "Pulsed" in
this design just means "thermostatically capped for safety," not
"low-duty-cycle," because the objective is now speed, not standby
efficiency.

If 1.7 kW/m is uncomfortably high for the actual hardware, the gap-vs-
power tradeoff table in `results_fast_melt/sleeve_sweep.csv` gives
cheaper options at a smaller gap (e.g. ~10 mm gap only needs ~1.0 kW/m,
~8 mm needs ~0.6 kW/m) — all still comfortably under the 3-minute
target with margin. Edit `REALISTIC_POWER_CAP`, `TARGET_S`,
`T_SAFETY_CAP` or the PCM dataset at the top of `fast_melt_design.py`
/ in `params.py` and rerun to fit an actual heater's power budget.

Output files (`results_fast_melt/`): `fig1_gap_sweeps.png` (the two
tradeoff curves), `fig2_liquid_fraction.png` (melt curve at the
recommended gap, with/without the sleeve), `fig3_temperature_profiles.png`
(radial snapshots — visibly two fronts, from the pipe and from the
sleeve, meeting in the middle), `fig4_sleeve_pulsing.png` (sleeve
temperature/power trace), `no_sleeve_sweep.csv`, `sleeve_sweep.csv`,
`summary.json`.

## Validation

This solver reuses, unmodified in form, the exact liquid-fraction and
secant-apparent-heat-capacity formulas from `PCM_MATLAB/liquidFraction.m`
and `PCM_MATLAB/apparentCapacity.m`, which are validated there against
the analytical 1-D Stefan problem to <0.5% melt-front error
(`PCM_MATLAB/validateStefan.m`). For this model specifically:
- **Grid convergence**: 120/240/400/600 radial cells agree on t99 to
  within ~1%.
- **Time-step convergence**: dt = 4/2/1/0.5 s agree on t99 to within
  ~0.1%.
- **Sanity check on the bisected minimum power**: a back-of-envelope
  energy balance (copper sensible heat + 1-D conduction penetration
  depth into the adjacent solid PCM over `t_ramp_target`) predicts the
  same order of magnitude (~200 W/m) as the numerically bisected
  value (291 W/m for the centered case).
