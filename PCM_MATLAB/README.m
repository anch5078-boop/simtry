%README  PCM_MATLAB -- Version 1 (conduction + phase change)
%
%   3D transient heat-conduction solver for a cylindrical PCM store
%   heated by two concentric copper cylindrical heater sleeves, with
%   uniform Cu-particle enhancement of the PCM and optional pulsed
%   heating. Conduction only (no natural convection) -- this is the
%   Version-1 model described in the project design notes; natural
%   convection is a planned higher-fidelity extension once this stage
%   is validated.
%
%   Runs unmodified in MATLAB or GNU Octave (tested on Octave 8.4).
%
%   This file is a plain MATLAB comment script, not a function --
%   open it in the editor to read it section by section, or run
%       publish('README.m')
%   to generate an HTML/PDF doc page from it.

%% Quick start
%
%   cd PCM_MATLAB
%   main             % builds geometry, runs the transient sim, plots + saves results/
%
%   Validate the solver against the analytical 1-D Stefan melting
%   problem:
%
%   validateStefan
%
%   Headless Octave note: if you run under octave-cli without a Qt
%   graphics toolkit, the fallback gnuplot backend does not honor
%   per-pixel image transparency, so the mid-height slice plots may
%   show the outside-the-cylinder region filled in rather than blank.
%   The underlying data is unaffected -- verified directly against
%   G.cylinderMask -- this is purely a headless-renderer limitation;
%   real MATLAB, or Octave with the qt toolkit, renders it correctly.
%
%   Run the heater-radius / particle-fraction case sweep (Section 13):
%
%   T = runParameterSweep([0.20 0.30 0.40], [0.55 0.70 0.85], [0 0.005 0.01 0.02]);
%   writetable(T, 'results/sweep.csv');

%% File structure
%
%   parameters.m           All physical/geometric/numerical constants
%                           (edit this to change a case)
%   createGeometry.m        Builds the Cartesian grid, cylindrical PCM
%                           mask, and the two heater-shell masks;
%                           reduces to an active-cell-only numbering
%                           for the solver
%   liquidFraction.m        Smoothed (cubic Hermite) liquid fraction
%                           f_l(T) and its derivative
%   apparentCapacity.m      Energy-conserving ("secant") apparent heat
%                           capacity -- see note below
%   compositeProperties.m   Maxwell effective conductivity
%                           (solid/liquid) + effective density/heat
%                           capacity for the Cu-particle-enhanced PCM
%   heaterPower.m           Pulsed/constant heater power -> volumetric
%                           heat generation Q''' in the heater cells
%   solveOneTimeStep.m      One implicit (backward-Euler + Picard)
%                           conduction/phase-change time step
%   runSimulation.m         Main time loop; tracks average liquid
%                           fraction, t50/t90/t99
%   plotResults.m           Standard result plots (liquid fraction
%                           history, temperature history, mid-height
%                           slices)
%   runParameterSweep.m     Batch runner over r1/R, r2/R, phi
%   validateStefan.m        1-D Stefan-problem validation (analytical
%                           vs. numerical melt front)
%   main.m                  Top-level driver script
%
%   Small shared utilities (kept as separate files rather than
%   script-local functions for MATLAB/Octave compatibility):
%   expandToFullGrid.m, frontPosition.m, formatTimeOrNA.m,
%   markMilestone.m, drawGuideline.m

%% Numerical method
%
%   Structured Cartesian finite-volume grid of cubic cells; only cells
%   with x^2+y^2 <= R^2 (and 0<=z<=H) are simulated ("active" cells)
%   -- the circular cross-section is staircase-approximated, which is
%   acceptable for a first screening model (quantify via a
%   grid-sensitivity study for the final paper). The outer cylindrical
%   boundary and domain edges are adiabatic by construction (no flux
%   term is ever added there).
%
%   The governing equation
%
%       C_app(T) dT/dt = div(k(T) grad T) + Q'''
%
%   is discretized implicitly (backward Euler) with harmonic-mean face
%   conductivities, and the resulting nonlinear system is handled with
%   a few Picard sub-iterations per step (P.picardIters, P.picardTol).

%% Energy-conserving phase-change capacity
%
%   The textbook apparent-heat-capacity model
%       C_app = rho*cp + rho*L*dfl/dT
%   uses the *local tangent* of the liquid-fraction curve. If a single
%   time step's temperature change is comparable to (or larger than)
%   the width of the smoothed mushy interval Tl-Ts, that tangent form
%   can "jump across" the latent-heat spike within one Picard iterate
%   and silently under-count the latent heat -- this was caught
%   directly by validateStefan.m (a naive implementation gave >600%
%   melt-front error). Instead, apparentCapacity.m uses the SECANT
%   slope of the enthalpy curve between the previous time level and
%   the current Picard iterate, which exactly conserves the
%   sensible+latent enthalpy change regardless of step size.
%   validateStefan.m confirms this: energy input through the boundary
%   matches stored enthalpy to <0.1% at any mesh/time-step resolution
%   tested.

%% Initializing exactly at Tm biases the T=Tm front diagnostic
%
%   One pitfall to be aware of if you extend validateStefan.m:
%   initializing the far-field/solid at exactly T=Tm places every
%   untouched cell precisely at the midpoint of the smoothed
%   liquid-fraction curve (f_l(Tm)=0.5 by construction) even though it
%   has absorbed no real heat. Comparing a "T crosses Tm" front
%   diagnostic against the sharp-interface analytical solution in that
%   state is biased by tens of percent (a first version of this
%   validation showed >600% error from this alone, before being traced
%   down) even though the scheme conserves energy exactly throughout.
%   validateStefan.m avoids this by initializing at Ts (the bottom of
%   the mushy interval, f_l=0 exactly) instead of Tm -- physically
%   indistinguishable from the one-phase assumption once dTmush is
%   small, and it removes the bias: the shipped script reproduces the
%   analytical melt front to <0.5% error (comfortably inside the
%   3-5% target band).

%% Cases (Section 12 of the design notes)
%
%   Case   PCM                          Heater
%   ----   --------------------------   ---------------------------------
%   C0     Pure PCM (phi=0)             single/base heater (set Ppeak on
%                                        heater 1 only, or r2i outside
%                                        the domain)
%   C1     Pure PCM (phi=0)             two cylindrical heaters
%   C2     Cu-particle PCM (phi>0)      two cylindrical heaters
%   C3     Cu-particle PCM (phi>0)      pulsed two-heater case
%                                        (P.pulseOn=true)
%
%   Set up each case by editing parameters.m (or overriding fields on
%   the struct returned by parameters() before calling createGeometry
%   / runSimulation), then compare t50, t90, t99, and R.TheaterMax(end)
%   across cases.
