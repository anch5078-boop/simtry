%PCM_MASTER_SIMULATION  Complete, single-file simulation of the "hot
%pipe / PCM / pulsed copper sleeve / PCM / insulated container wall"
%problem: 1-D radial (cylindrical) transient conduction + phase-change
%(enthalpy method), covering both:
%
%   Part 1 -- the geometry as drawn (1" pipe, 32.3 mm pipe->wall gap,
%             copper sleeve centered in the gap, pulse-heated at the
%             lowest power that reaches the 40 C PCM melt point within
%             300 s), plus a no-sleeve reference and a sleeve-position
%             sweep to find the minimum achievable melting time for
%             this architecture (~2h20m, sleeve at ~85% of the gap).
%
%   Part 2 -- the same architecture re-dimensioned so the whole PCM
%             charge melts in under 3 minutes (pipe radius and PCM
%             material held fixed): a 12 mm gap, sleeve at 70% of the
%             gap, 1700 W/m, safety-capped at 180 C -- melts in ~167 s.
%
%THIS SINGLE FILE contains everything: every helper function this
%model needs (geometry, the smoothed-liquid-fraction / energy-
%conserving secant apparent-heat-capacity phase-change scheme,
%backward-Euler + Picard implicit solver, bisection helpers, plotting)
%followed by the driver code that actually runs both parts. Nothing
%else is required.
%
%   Run directly (MATLAB, or Octave via `octave-cli PCM_master_simulation.m`):
%       >> PCM_master_simulation
%
%   Total run time: a few minutes (most of it Part 2's gap-width x
%   power bisection sweep). Results (plots, CSVs, summary.txt) are
%   written to ./results_master/asDrawn/ and ./results_master/fastMelt/.
%
%   To change the case: edit parametersRadial() below (Part 1's
%   geometry/material/control constants) or the TARGET_S /
%   T_SAFETY_CAP / REALISTIC_POWER_CAP constants at the top of the
%   Part 2 driver section further down. Nothing else needs to change.
%
%   Numerical method: smoothed cubic-Hermite liquid fraction; energy-
%   conserving "secant" apparent heat capacity (exactly conserves
%   sensible+latent enthalpy regardless of time-step size -- the naive
%   tangent-derivative form can under-count latent heat and is NOT
%   used here); harmonic-mean face conductivity (correct for the
%   PCM/copper interface where k differs ~2000x); implicit backward-
%   Euler time stepping with Picard sub-iterations for the k(T)/Capp(T)
%   nonlinearity. Dirichlet inner boundary (pipe surface, fixed at the
%   hot-water temperature -- pipe-wall resistance neglected); adiabatic
%   outer boundary (insulated container wall). See this repo's
%   PCM_MATLAB_radial/README.md and PCM_python/README.md for the full
%   write-up, assumptions table, and validation checks (this file's
%   logic is cross-validated against an independent Python
%   implementation of the same discretization -- results match to the
%   last digit).
%
%   Tested on GNU Octave 8.4. Octave note: local functions in a script
%   file must be defined before their first use in the file, hence the
%   `1;` marker below and every function appearing before the driver
%   code that calls them (MATLAB allows either order; this file's order
%   works in both).

1;


% ------------------------------------------------------------------------
% source: parametersRadial.m
% ------------------------------------------------------------------------
function P = parametersRadial()
%PARAMETERSRADIAL  Every physical, geometric and control constant for the
%"hot pipe / PCM / pulsed copper sleeve / PCM / insulated wall" case, as
%drawn (sleeve centered in the pipe->wall gap). This is the MATLAB/Octave
%counterpart of PCM_python/params.py -- see that file's comments for the
%reasoning behind each assumed value; nothing here should need to change
%except by editing this file.

%% ---- Geometry (metres) -- rPipe is given, the rest is ASSUMED --------
P.rPipe   = 0.0127;    % 1" pipe OD / PCM-facing radius, given         [m]
P.rWall   = 0.045;     % container inner radius (ASSUMED)              [m]
P.tSleeve = 0.0015;    % copper sleeve radial thickness (ASSUMED)      [m]
P.centerFrac = 0.5;    % sleeve mid-radius as a fraction of the pipe->
                       % wall gap; 0.5 = "centered in the gap" as drawn

%% ---- Grid -------------------------------------------------------------
P.nCells = 240;

%% ---- PCM properties (paraffin/RT42-class wax, Tm ~ 40 C) -- ASSUMED,
%% replace with your PCM's real datasheet values.
P.rhoPCM = 830;        % kg/m^3 (single value, both phases)
P.cpPCM  = 2100;       % J/(kg K)
P.kPCMs  = 0.20;       % W/(m K), solid
P.kPCMl  = 0.15;       % W/(m K), liquid
P.LPCM   = 170000;     % J/kg, latent heat of fusion
P.Tm     = 40;         % deg C, nominal melt point -- given
P.dTmush = 4;          % deg C, smoothed mushy-zone width (numerical
                       % regularization only)

%% ---- Copper (sleeve) ----------------------------------------------------
P.rhoCu = 8960;
P.cpCu  = 385;
P.kCu   = 400;

%% ---- Boundary / operating conditions -- given --------------------------
P.Thotwater = 120;     % deg C, Dirichlet at the pipe surface (pipe wall
                       % resistance neglected)
P.Tinit     = 20;      % deg C, uniform initial temperature (ASSUMED)

%% ---- Sleeve pulse-heating control ---------------------------------------
% "Pulse heated at lowest power to reach 40 C": bisect for the smallest
% sleeve peak power that reaches Tm within tRampTarget seconds of a cold
% start, then hold there with a bang-bang thermostat (see
% runFullMeltRadial.m).
P.tRampTarget = 300;   % s (ASSUMED design target)

%% ---- Time stepping -------------------------------------------------------
P.dt = 2;              % s
P.tMax = 6*3600;       % s, safety cap
P.picardIters = 4;
P.picardTol   = 1e-3;

end


% ------------------------------------------------------------------------
% source: materialStruct.m
% ------------------------------------------------------------------------
function matp = materialStruct(P)
%MATERIALSTRUCT  Pack the scalar material properties solveOneStepRadial
%needs (mushy-zone bounds + effective rho*cp / rho*L) from the parameter
%struct P. Kept as a tiny separate step so P itself stays the single
%source of truth (mirrors PCM_python/params.py's Params.material_dict()).

matp.Ts = P.Tm - P.dTmush/2;
matp.Tl = P.Tm + P.dTmush/2;
matp.kPCMs = P.kPCMs;
matp.kPCMl = P.kPCMl;
matp.rhocpPCM = P.rhoPCM * P.cpPCM;
matp.rhoLPCM  = P.rhoPCM * P.LPCM;
matp.kCu = P.kCu;
matp.rhocpCu = P.rhoCu * P.cpCu;
matp.Thotwater = P.Thotwater;

end


% ------------------------------------------------------------------------
% source: sleeveRadii.m
% ------------------------------------------------------------------------
function [rIn, rOut] = sleeveRadii(P, centerFrac)
%SLEEVERADII  Inner/outer radius of the copper sleeve for a given
%fractional position across the pipe->wall gap (0 = touching the pipe,
%1 = touching the wall). Defaults to P.centerFrac ("centered in the
%gap", as drawn) when centerFrac is omitted.

if nargin < 2 || isempty(centerFrac)
    centerFrac = P.centerFrac;
end

rMid = P.rPipe + centerFrac * (P.rWall - P.rPipe);
rIn  = rMid - P.tSleeve/2;
rOut = rMid + P.tSleeve/2;

end


% ------------------------------------------------------------------------
% source: createGeometryRadial.m
% ------------------------------------------------------------------------
function G = createGeometryRadial(P, centerFrac)
%CREATEGEOMETRYRADIAL  Build the uniform-spacing 1-D radial finite-volume
%grid from P.rPipe to P.rWall, and the copper-sleeve cell mask.
%
%   G = CREATEGEOMETRYRADIAL(P, centerFrac)
%
%       P          : parameter struct from parametersRadial.m
%       centerFrac : (optional) sleeve mid-radius as a fraction of the
%                    pipe->wall gap; defaults to P.centerFrac.
%
%   All of G's "volume"/"area" fields are per unit AXIAL length of pipe
%   (the standard reduction of an axisymmetric, axially-uniform problem
%   to one radial dimension): V(i) is really an area [m^2], Aface(i) is
%   really a circumference [m].

if nargin < 2
    centerFrac = P.centerFrac;
end

N = P.nCells;
rFace = linspace(P.rPipe, P.rWall, N+1)';
rc    = 0.5*(rFace(1:end-1) + rFace(2:end));
dr    = rFace(2) - rFace(1);
V     = pi*(rFace(2:end).^2 - rFace(1:end-1).^2);
Aface = 2*pi*rFace;

[rSleeveIn, rSleeveOut] = sleeveRadii(P, centerFrac);
isCu = rc >= rSleeveIn & rc <= rSleeveOut;

if ~any(isCu)
    error('createGeometryRadial:noSleeve', ...
        'No grid cell falls inside the sleeve band [%.5g, %.5g] m -- increase nCells or tSleeve.', ...
        rSleeveIn, rSleeveOut);
end

G.N = N;
G.rFace = rFace;
G.rc = rc;
G.dr = dr;
G.V = V;
G.Aface = Aface;
G.isCu = isCu;
G.isPcm = ~isCu;
G.rSleeveIn = rSleeveIn;
G.rSleeveOut = rSleeveOut;

end


% ------------------------------------------------------------------------
% source: createGeometryNoSleeve.m
% ------------------------------------------------------------------------
function G = createGeometryNoSleeve(P)
%CREATEGEOMETRYNOSLEEVE  Same radial grid as createGeometryRadial.m but
%with no copper sleeve anywhere (pure PCM from pipe to wall) -- the
%"pipe conduction only" reference case.

N = P.nCells;
rFace = linspace(P.rPipe, P.rWall, N+1)';
rc    = 0.5*(rFace(1:end-1) + rFace(2:end));
dr    = rFace(2) - rFace(1);
V     = pi*(rFace(2:end).^2 - rFace(1:end-1).^2);
Aface = 2*pi*rFace;

G.N = N;
G.rFace = rFace;
G.rc = rc;
G.dr = dr;
G.V = V;
G.Aface = Aface;
G.isCu  = false(N,1);
G.isPcm = true(N,1);
G.rSleeveIn = NaN;
G.rSleeveOut = NaN;

end


% ------------------------------------------------------------------------
% source: liquidFractionRadial.m
% ------------------------------------------------------------------------
function fl = liquidFractionRadial(T, Ts, Tl)
%LIQUIDFRACTIONRADIAL  Smooth (cubic Hermite / "smoothstep") liquid
%fraction, identical in form to PCM_MATLAB/liquidFraction.m:
%
%       fl(T) = 0                          , T <= Ts
%              = 3*xi^2 - 2*xi^3            , Ts < T < Tl
%              = 1                          , T >= Tl
%   where xi = (T - Ts) / (Tl - Ts).
%
%   T may be any array; fl is returned the same size/shape.

dT = Tl - Ts;
if dT <= 0
    error('liquidFractionRadial:badInterval', 'Require Tl > Ts.');
end

xi = (T - Ts) / dT;
xi = min(max(xi, 0), 1);

fl = 3*xi.^2 - 2*xi.^3;

end


% ------------------------------------------------------------------------
% source: liquidFractionDerivRadial.m
% ------------------------------------------------------------------------
function dfldT = liquidFractionDerivRadial(T, Ts, Tl)
%LIQUIDFRACTIONDERIVRADIAL  d(fl)/dT for the smoothed liquid fraction in
%liquidFractionRadial.m -- used only as the tiny-step fallback inside
%apparentCapacityRadial.m.

dT = Tl - Ts;
xi = min(max((T - Ts) / dT, 0), 1);
dfldT = 6*xi.*(1 - xi) / dT;

outside = (T <= Ts) | (T >= Tl);
dfldT(outside) = 0;

end


% ------------------------------------------------------------------------
% source: apparentCapacityRadial.m
% ------------------------------------------------------------------------
function Capp = apparentCapacityRadial(Tguess, Tn, Ts, Tl, rhocp, rhoL)
%APPARENTCAPACITYRADIAL  Energy-conserving ("secant") apparent heat
%capacity -- identical in form to PCM_MATLAB/apparentCapacity.m (see
%that file's header comment for the full derivation and for why the
%naive tangent-derivative form under-counts latent heat when a step's
%dT is comparable to the mushy-zone width):
%
%   Capp = rhocp + rhoL * (fl(Tguess)-fl(Tn)) / (Tguess-Tn)
%
%falling back to the tangent derivative dfl/dT where Tguess and Tn
%coincide to machine precision.

dT = Tguess - Tn;
tinyStep = abs(dT) < 1e-8;

Capp = zeros(size(Tguess));

if any(~tinyStep(:))
    flG = liquidFractionRadial(Tguess(~tinyStep), Ts, Tl);
    flN = liquidFractionRadial(Tn(~tinyStep),     Ts, Tl);
    Capp(~tinyStep) = rhocp + rhoL .* (flG - flN) ./ dT(~tinyStep);
end

if any(tinyStep(:))
    dfldT = liquidFractionDerivRadial(Tguess(tinyStep), Ts, Tl);
    Capp(tinyStep) = rhocp + rhoL .* dfldT;
end

end


% ------------------------------------------------------------------------
% source: solveOneStepRadial.m
% ------------------------------------------------------------------------
function [Tnew, fl] = solveOneStepRadial(Tn, G, matp, Qtotal, dt, picardIters, picardTol)
%SOLVEONESTEPRADIAL  Advance the radial temperature field by one implicit
%(backward-Euler) time step on the grid G, solving
%
%   Capp(T) * V/dt * (T - Tn) = sum_faces G_face*(T_nb - T) + Q_cell
%
%with a Dirichlet inner boundary (fixed at matp.Thotwater, the pipe
%surface) and an adiabatic outer boundary (the insulated container
%wall, reproduced simply by never adding a flux term there). k(T) and
%Capp(T) are nonlinear through the liquid fraction, so a few Picard
%sub-iterations are used, exactly as in PCM_MATLAB/solveOneTimeStep.m.
%
%   [Tnew, fl] = SOLVEONESTEPRADIAL(Tn, G, matp, Qtotal, dt, picardIters, picardTol)
%
%       Tn      : N x 1 temperature at the start of the step   [deg C]
%       G       : geometry struct from createGeometryRadial.m / createGeometryNoSleeve.m
%       matp    : material struct from materialStruct.m
%       Qtotal  : instantaneous sleeve power this step          [W/m of pipe length]
%       dt      : time step                                     [s]
%
%       Tnew, fl : N x 1 temperature / liquid fraction at the end of the step

N = G.N;
V = G.V;
dr = G.dr;
isCu  = G.isCu;
isPcm = G.isPcm;

VCu = sum(V(isCu));
qCell = zeros(N, 1);
if VCu > 0 && Qtotal ~= 0
    qCell(isCu) = Qtotal * V(isCu) / VCu;   % distribute by volume
end

Tguess = Tn;
idx = (1:N)';

for it = 1:picardIters
    fl = liquidFractionRadial(Tguess, matp.Ts, matp.Tl);

    k = zeros(N, 1);
    k(isPcm) = (1 - fl(isPcm)).*matp.kPCMs + fl(isPcm).*matp.kPCMl;
    k(isCu)  = matp.kCu;

    Capp = zeros(N, 1);
    Capp(isPcm) = apparentCapacityRadial(Tguess(isPcm), Tn(isPcm), ...
                                          matp.Ts, matp.Tl, matp.rhocpPCM, matp.rhoLPCM);
    Capp(isCu)  = matp.rhocpCu;

    % Harmonic-mean interior face conductance (uniform dr -> cell
    % centres equidistant from the shared face).
    kFace = 2*k(1:end-1).*k(2:end) ./ (k(1:end-1) + k(2:end));
    Gint  = kFace .* G.Aface(2:end-1) / dr;         % N-1 x 1

    % Inner Dirichlet boundary conductance (pipe surface, dr/2 from the
    % first cell centre).
    GbcIn = 2*k(1)*G.Aface(1) / dr;

    diagMain = Capp .* V / dt;
    diagMain(1:end-1) = diagMain(1:end-1) + Gint;
    diagMain(2:end)   = diagMain(2:end)   + Gint;
    diagMain(1)       = diagMain(1) + GbcIn;

    offDiag = -Gint;

    rows = [idx(1:end-1); idx(2:end); idx];
    cols = [idx(2:end); idx(1:end-1); idx];
    vals = [offDiag; offDiag; diagMain];
    A = sparse(rows, cols, vals, N, N);

    b = Capp .* V / dt .* Tn + qCell;
    b(1) = b(1) + GbcIn * matp.Thotwater;

    Tnew = A \ b;

    change = max(abs(Tnew - Tguess));
    scale  = max(1, max(abs(Tguess)));
    Tguess = Tnew;

    if change < picardTol * scale
        break;
    end
end

fl = liquidFractionRadial(Tnew, matp.Ts, matp.Tl);

end


% ------------------------------------------------------------------------
% source: volumeAvg.m
% ------------------------------------------------------------------------
function avg = volumeAvg(field, V, mask)
%VOLUMEAVG  Volume-weighted average of FIELD over cells where MASK is
%true (cell volumes -- really per-unit-length areas -- vary with radius
%in this cylindrical grid, so a plain mean would be wrong).

avg = sum(field(mask).*V(mask)) / sum(V(mask));

end


% ------------------------------------------------------------------------
% source: rampFinalSleeveTemp.m
% ------------------------------------------------------------------------
function Tend = rampFinalSleeveTemp(P, G, matp, Ppeak, tRamp)
%RAMPFINALSLEEVETEMP  Run the sleeve-heating ramp (heater full-on at
%Ppeak from a cold start, pipe boundary already at Thotwater throughout)
%for tRamp seconds and return the volume-averaged sleeve temperature at
%the end. Used as the scalar objective for the minimum-power bisection
%in findMinPowerRadial.m.

N = G.N;
T = P.Tinit * ones(N, 1);
nSteps = max(1, round(tRamp / P.dt));

for i = 1:nSteps
    T = solveOneStepRadial(T, G, matp, Ppeak, P.dt, P.picardIters, P.picardTol);
end

Tend = volumeAvg(T, G.V, G.isCu);

end


% ------------------------------------------------------------------------
% source: findMinPowerRadial.m
% ------------------------------------------------------------------------
function [Pmin, Tend] = findMinPowerRadial(P, G, matp, tolTemp, maxIter)
%FINDMINPOWERRADIAL  Bisect for the smallest sleeve peak power [W/m of
%pipe length] such that the volume-averaged sleeve temperature reaches
%P.Tm within P.tRampTarget seconds of a cold start (pipe boundary
%already at Thotwater throughout, as in the full problem). This is the
%literal "pulse heated at lowest power to reach 40 C" rule.

if nargin < 4 || isempty(tolTemp), tolTemp = 0.05; end
if nargin < 5 || isempty(maxIter), maxIter = 40; end

pLo = 0.0;
Tend = rampFinalSleeveTemp(P, G, matp, pLo, P.tRampTarget);
if Tend >= P.Tm
    Pmin = pLo;   % even zero power gets there within the ramp window
    return;
end

pHi = 50.0;
found = false;
for i = 1:30
    Tend = rampFinalSleeveTemp(P, G, matp, pHi, P.tRampTarget);
    if Tend >= P.Tm
        found = true;
        break;
    end
    pHi = pHi * 2.0;
end
if ~found
    error('findMinPowerRadial:noBracket', 'Could not bracket a sufficient sleeve power.');
end

lo = pLo; hi = pHi;
for i = 1:maxIter
    mid = 0.5*(lo+hi);
    Tend = rampFinalSleeveTemp(P, G, matp, mid, P.tRampTarget);
    if Tend >= P.Tm
        hi = mid;
    else
        lo = mid;
    end
    if abs(Tend - P.Tm) < tolTemp
        break;
    end
end
Pmin = hi;
Tend = rampFinalSleeveTemp(P, G, matp, Pmin, P.tRampTarget);

end


% ------------------------------------------------------------------------
% source: runFullMeltRadial.m
% ------------------------------------------------------------------------
function R = runFullMeltRadial(P, G, matp, Ppeak, Tsetpoint, recordEvery, storeSnapshots)
%RUNFULLMELTRADIAL  Full transient melt simulation: bang-bang thermostat
%holds the sleeve at Tsetpoint (heater ON at Ppeak while the volume-avg
%sleeve temperature is below Tsetpoint, OFF otherwise) while the pipe
%boundary sits fixed at Thotwater throughout. Tsetpoint defaults to
%P.Tm (the efficiency-oriented "lowest power to just reach/hold the
%melt point" rule); pass a higher value for a speed-oriented "run the
%sleeve as hard as the supply allows, up to a safety cutoff" case.
%
%   R = RUNFULLMELTRADIAL(P, G, matp, Ppeak, Tsetpoint, recordEvery, storeSnapshots)
%
%   Returns a struct R with time histories (R.t, R.favg, R.Tsleeve,
%   R.Tmax, R.power, R.dutyOn), the t50/t90/t99/t999 (99.9%) melt
%   milestones based on the volume-weighted average PCM liquid
%   fraction, the overall sleeve duty cycle, and (if storeSnapshots)
%   full T(r)/fl(r) snapshots at each recorded step.

if nargin < 5 || isempty(Tsetpoint),     Tsetpoint = P.Tm; end
if nargin < 6 || isempty(recordEvery),   recordEvery = 30; end
if nargin < 7 || isempty(storeSnapshots), storeSnapshots = false; end

N = G.N;
V = G.V;
isCu  = G.isCu;
isPcm = G.isPcm;
hasSleeve = any(isCu);

T = P.Tinit * ones(N, 1);
t = 0;
nSteps = round(P.tMax / P.dt);

tHist     = zeros(nSteps+1, 1);
favgHist  = zeros(nSteps+1, 1);
TsHist    = zeros(nSteps+1, 1);
TmaxHist  = zeros(nSteps+1, 1);
powHist   = zeros(nSteps+1, 1);
dutyHist  = false(nSteps+1, 1);

tHist(1) = 0; favgHist(1) = 0;
TsHist(1) = P.Tinit; TmaxHist(1) = P.Tinit;

snapshots = struct('t', {}, 'T', {}, 'fl', {});
if storeSnapshots
    snapshots(end+1) = struct('t', 0, 'T', T, 'fl', zeros(N,1));
end

t50 = NaN; t90 = NaN; t99 = NaN; t999 = NaN;
heaterOnTime = 0;
nRec = 1;

for step = 1:nSteps
    if hasSleeve
        TsleeveNow = volumeAvg(T, V, isCu);
        on = TsleeveNow < Tsetpoint;
    else
        on = false;
    end
    Papplied = 0;
    if on
        Papplied = Ppeak;
        heaterOnTime = heaterOnTime + P.dt;
    end

    [T, fl] = solveOneStepRadial(T, G, matp, Papplied, P.dt, P.picardIters, P.picardTol);
    t = t + P.dt;

    favg = volumeAvg(fl, V, isPcm);

    if mod(step, recordEvery) == 0
        nRec = nRec + 1;
        tHist(nRec) = t;
        favgHist(nRec) = favg;
        if hasSleeve
            TsHist(nRec) = volumeAvg(T, V, isCu);
        else
            TsHist(nRec) = NaN;
        end
        TmaxHist(nRec) = max(T);
        powHist(nRec) = Papplied;
        dutyHist(nRec) = on;
        if storeSnapshots
            snapshots(end+1) = struct('t', t, 'T', T, 'fl', fl); %#ok<AGROW>
        end
    end

    if isnan(t50)  && favg >= 0.50,  t50  = t; end
    if isnan(t90)  && favg >= 0.90,  t90  = t; end
    if isnan(t99)  && favg >= 0.99,  t99  = t; end
    if isnan(t999) && favg >= 0.999, t999 = t; end

    if ~isnan(t999)
        break;
    end
end

tHist    = tHist(1:nRec);
favgHist = favgHist(1:nRec);
TsHist   = TsHist(1:nRec);
TmaxHist = TmaxHist(1:nRec);
powHist  = powHist(1:nRec);
dutyHist = dutyHist(1:nRec);

R.t = tHist; R.favg = favgHist; R.Tsleeve = TsHist; R.Tmax = TmaxHist;
R.power = powHist; R.dutyOn = dutyHist;
R.t50 = t50; R.t90 = t90; R.t99 = t99; R.t999 = t999;
R.tFinal = t;
if t > 0
    R.dutyCycle = heaterOnTime / t;
else
    R.dutyCycle = NaN;
end
R.converged = ~isnan(t999);
R.favgFinal = favgHist(end);
R.snapshots = snapshots;

end


% ------------------------------------------------------------------------
% source: runPositionSweepRadial.m
% ------------------------------------------------------------------------
function rows = runPositionSweepRadial(P, matp, fracs, recordEvery, verbose)
%RUNPOSITIONSWEEPRADIAL  Vary the copper sleeve's radial position across
%the pipe->wall gap (holding the "lowest power to reach Tm within the
%ramp target" control rule fixed at each position) to find the
%placement that minimizes total melting time.
%
%   rows = RUNPOSITIONSWEEPRADIAL(P, matp, fracs, recordEvery, verbose)
%
%   Returns an array of structs, one per fraction tried, each with
%   frac, rIn, rOut, rMid, Pmin, t50, t90, t99, t999, dutyCycle.

if nargin < 3 || isempty(fracs), fracs = linspace(0.15, 0.85, 8); end
if nargin < 4 || isempty(recordEvery), recordEvery = 60; end
if nargin < 5 || isempty(verbose), verbose = true; end

rows = struct('frac', {}, 'rIn', {}, 'rOut', {}, 'rMid', {}, 'Pmin', {}, ...
              't50', {}, 't90', {}, 't99', {}, 't999', {}, 'dutyCycle', {}, 'converged', {});

for i = 1:numel(fracs)
    frac = fracs(i);
    G = createGeometryRadial(P, frac);
    [Pmin, ~] = findMinPowerRadial(P, G, matp);
    R = runFullMeltRadial(P, G, matp, Pmin, P.Tm, recordEvery, false);

    row.frac = frac;
    row.rIn = G.rSleeveIn; row.rOut = G.rSleeveOut;
    row.rMid = 0.5*(G.rSleeveIn + G.rSleeveOut);
    row.Pmin = Pmin;
    row.t50 = R.t50; row.t90 = R.t90; row.t99 = R.t99; row.t999 = R.t999;
    row.dutyCycle = R.dutyCycle;
    row.converged = R.converged;
    rows(end+1) = row; %#ok<AGROW>

    if verbose
        fprintf(['  frac=%4.2f  rMid=%6.2fmm  Pmin=%7.2f W/m  ', ...
                 't50=%8s  t90=%8s  t99=%8s  t999=%8s\n'], ...
            frac, row.rMid*1000, Pmin, ...
            numOrNA(R.t50), numOrNA(R.t90), numOrNA(R.t99), numOrNA(R.t999));
    end
end

end

function s = numOrNA(x)
if isnan(x)
    s = 'n/a';
else
    s = sprintf('%.1f', x);
end
end


% ------------------------------------------------------------------------
% source: t999ForPower.m
% ------------------------------------------------------------------------
function [t999, R] = t999ForPower(P, G, matp, power, Tsetpoint, recordEvery)
%T999FORPOWER  Convenience wrapper: run runFullMeltRadial at a given
%sleeve peak power (bang-bang up to Tsetpoint) and return just the
%full-melt (99.9%) milestone, plus the full result struct.

if nargin < 6 || isempty(recordEvery), recordEvery = 100; end
R = runFullMeltRadial(P, G, matp, power, Tsetpoint, recordEvery, false);
t999 = R.t999;

end


% ------------------------------------------------------------------------
% source: minPowerForTarget.m
% ------------------------------------------------------------------------
function [Pmin, t999] = minPowerForTarget(P, G, matp, targetS, Tsetpoint, pHi0, pCap, maxIter)
%MINPOWERFORTARGET  Bisect for the smallest sleeve peak power such that
%the full-melt time (t999) is <= targetS, with the sleeve bang-bang
%thermostat capped at Tsetpoint (a safety ceiling, not an efficiency
%target -- see runFullMeltRadial.m). Returns [] (NaN) if no power up to
%pCap achieves the target.

if nargin < 6 || isempty(pHi0),  pHi0 = 50.0;    end
if nargin < 7 || isempty(pCap),  pCap = 20000.0; end
if nargin < 8 || isempty(maxIter), maxIter = 28; end

pHi = pHi0;
while true
    t999 = t999ForPower(P, G, matp, pHi, Tsetpoint);
    if ~isnan(t999) && t999 <= targetS
        break;
    end
    pHi = pHi * 1.6;
    if pHi > pCap
        Pmin = NaN; t999 = NaN;
        return;
    end
end

pLo = 0.0;
for i = 1:maxIter
    mid = 0.5*(pLo + pHi);
    t999 = t999ForPower(P, G, matp, mid, Tsetpoint);
    if ~isnan(t999) && t999 <= targetS
        pHi = mid;
    else
        pLo = mid;
    end
end
Pmin = pHi;
t999 = t999ForPower(P, G, matp, Pmin, Tsetpoint);

end


% ------------------------------------------------------------------------
% source: positionSweepQuick.m
% ------------------------------------------------------------------------
function bestFrac = positionSweepQuick(P, matp, fracs, refPower, Tsetpoint)
%POSITIONSWEEPQUICK  Rank sleeve positions by full-melt time at a fixed
%reference power (just to locate the local optimum quickly, before the
%real minimum-power bisection is run at that position) -- mirrors
%PCM_python/fast_melt_design.py's position_sweep().

if nargin < 3 || isempty(fracs), fracs = linspace(0.3, 0.85, 12); end
if nargin < 4 || isempty(refPower), refPower = 1500.0; end
if nargin < 5 || isempty(Tsetpoint), Tsetpoint = 180.0; end

bestT999 = Inf;
bestFrac = fracs(1);
for i = 1:numel(fracs)
    G = createGeometryRadial(P, fracs(i));
    t999 = t999ForPower(P, G, matp, refPower, Tsetpoint);
    if ~isnan(t999) && t999 < bestT999
        bestT999 = t999;
        bestFrac = fracs(i);
    end
end

end


% ------------------------------------------------------------------------
% source: sweepNoSleeveGaps.m
% ------------------------------------------------------------------------
function rows = sweepNoSleeveGaps(rPipe, gapListMm, nCells, dt, tMax, Pbase, verbose)
%SWEEPNOSLEEVEGAPS  Pure pipe-conduction (no sleeve) full-melt time vs.
%PCM annulus gap width -- shows the "no auxiliary heater needed" limit
%for the <3-minute target.

if nargin < 7 || isempty(verbose), verbose = true; end

rows = struct('gapMm', {}, 't999', {});
for i = 1:numel(gapListMm)
    gapMm = gapListMm(i);
    P = Pbase;
    P.rPipe = rPipe;
    P.rWall = rPipe + gapMm/1000;
    P.nCells = nCells; P.dt = dt; P.tMax = tMax;
    matp = materialStruct(P);
    G = createGeometryNoSleeve(P);
    R = runFullMeltRadial(P, G, matp, 0, P.Tm, 40, false);

    row.gapMm = gapMm; row.t999 = R.t999;
    rows(end+1) = row; %#ok<AGROW>
    if verbose
        fprintf('  [no sleeve] gap=%3d mm  t999=%s\n', gapMm, numOrNAstr(R.t999));
    end
end

end

function s = numOrNAstr(x)
if isnan(x), s = 'n/a'; else, s = sprintf('%.2f', x); end
end


% ------------------------------------------------------------------------
% source: sweepWithSleeveGaps.m
% ------------------------------------------------------------------------
function rows = sweepWithSleeveGaps(rPipe, gapListMm, nCells, dt, tMax, tSleeve, Pbase, targetS, Tsetpoint, verbose)
%SWEEPWITHSLEEVEGAPS  With the pulsed copper sleeve: for each PCM
%annulus gap width, re-optimize the sleeve's radial position (quick
%ranking at a fixed reference power) and bisect for the minimum sleeve
%peak power that still melts everything within targetS seconds
%(sleeve bang-bang capped at Tsetpoint for safety).

if nargin < 8  || isempty(targetS),   targetS = 170.0; end
if nargin < 9  || isempty(Tsetpoint), Tsetpoint = 180.0; end
if nargin < 10 || isempty(verbose),   verbose = true; end

rows = struct('gapMm', {}, 'frac', {}, 'Pmin', {}, 't999', {});
for i = 1:numel(gapListMm)
    gapMm = gapListMm(i);
    P = Pbase;
    P.rPipe = rPipe;
    P.rWall = rPipe + gapMm/1000;
    P.nCells = nCells; P.dt = dt; P.tMax = tMax; P.tSleeve = tSleeve;
    matp = materialStruct(P);

    frac = positionSweepQuick(P, matp, linspace(0.3, 0.85, 12), 1500.0, Tsetpoint);
    G = createGeometryRadial(P, frac);
    [Pmin, t999] = minPowerForTarget(P, G, matp, targetS, Tsetpoint);

    row.gapMm = gapMm; row.frac = frac; row.Pmin = Pmin; row.t999 = t999;
    rows(end+1) = row; %#ok<AGROW>
    if verbose
        fprintf('  [with sleeve] gap=%3d mm  frac=%.2f  Pmin=%s  t999=%s\n', ...
            gapMm, frac, numOrNAstr(Pmin), numOrNAstr(t999));
    end
end

end



% ------------------------------------------------------------------------
% source: numOrEmpty.m
% ------------------------------------------------------------------------
function s = numOrEmpty(x)
%NUMOREMPTY  "%.4f"-formatted number, or an empty string for NaN
%(kept as a separate file, not a script-local function, for
%MATLAB/Octave compatibility -- mirrors the convention already used by
%PCM_MATLAB/formatTimeOrNA.m et al.).
if isnan(x)
    s = '';
else
    s = sprintf('%.4f', x);
end
end


% ------------------------------------------------------------------------
% source: formatTimeOrNARadial.m
% ------------------------------------------------------------------------
function s = formatTimeOrNARadial(x)
%FORMATTIMEORNARADIAL  "Hh Mm Ss (N s)" for a finite scalar, else "n/a".
if isempty(x) || isnan(x)
    s = 'n/a';
    return;
end
h = floor(x/3600);
m = floor(mod(x,3600)/60);
sec = mod(x,60);
s = sprintf('%dh %02dm %02ds (%.0f s)', h, m, sec, x);
end


% ------------------------------------------------------------------------
% source: makePlotsAsDrawn.m
% ------------------------------------------------------------------------
function makePlotsAsDrawn(P, Rbase, Rref, Ropt, rowsSorted, best, bandBase, bandOpt, outDir)
%MAKEPLOTSASDRAWN  Result plots for mainAsDrawn.m -- mirrors
%PCM_python/main.py's make_plots(): melting-progress comparison,
%position-sweep curve, temperature-profile snapshots, sleeve pulsing
%trace.

% ---- fig1: liquid fraction history -----------------------------------
f = figure('Visible', 'off');
plot(Rbase.t/60, Rbase.favg*100, '-', 'LineWidth', 1.8); hold on;
plot(Ropt.t/60,  Ropt.favg*100,  '-', 'LineWidth', 1.8);
plot(Rref.t/60,  Rref.favg*100,  '--', 'LineWidth', 1.8);
yline_manual(99);
xlabel('time [min]'); ylabel('PCM average liquid fraction [%]');
title('Melting progress: effect of the pulsed copper sleeve');
legend('Sleeve centered in gap (as drawn)', ...
       sprintf('Optimal position (%.0f%% pipe->wall)', best.frac*100), ...
       'No sleeve (pipe conduction only)', 'Location', 'southeast');
grid on;
print(f, fullfile(outDir, 'fig1_liquid_fraction.png'), '-dpng', '-r150');
close(f);

% ---- fig2: sweep curve -------------------------------------------------
fracs = [rowsSorted.frac];
t99s  = [rowsSorted.t99]/60;
t999s = [rowsSorted.t999]/60;
f = figure('Visible', 'off');
plot(fracs, t99s, 'o-', 'MarkerSize', 4); hold on;
plot(fracs, t999s, 's-', 'MarkerSize', 4);
xline_manual(best.frac, 'g');
xline_manual(0.5, [0.5 0.5 0.5]);
xlabel('sleeve position: fraction of pipe->wall gap');
ylabel('melting time [min]');
title('Melting time vs. copper-sleeve radial position');
legend('t99 (99% melted)', 't99.9 (fully melted)', 'optimum', 'centered (as drawn)', ...
       'Location', 'northeast');
grid on;
print(f, fullfile(outDir, 'fig2_position_sweep.png'), '-dpng', '-r150');
close(f);

% ---- fig3: temperature profile snapshots -------------------------------
G = createGeometryRadial(P, best.frac);
snaps = Ropt.snapshots;
nSnap = numel(snaps);
pick = round(linspace(1, nSnap, min(7, nSnap)));
f = figure('Visible', 'off'); hold on;
cmap = jet(numel(pick));
legendEntries = {};
for i = 1:numel(pick)
    s = snaps(pick(i));
    plot(G.rc*1000, s.T, 'Color', cmap(i,:), 'LineWidth', 1.4);
    legendEntries{end+1} = sprintf('t=%.0f min', s.t/60); %#ok<AGROW>
end
plot(xlim(), [P.Tm P.Tm], 'r--', 'LineWidth', 1);
legendEntries{end+1} = 'PCM melt point';
patch([bandOpt(1) bandOpt(2) bandOpt(2) bandOpt(1)]*1000, ...
      [min(ylim()) min(ylim()) max(ylim()) max(ylim())], ...
      [1 0.65 0], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
legendEntries{end+1} = 'copper sleeve';
xlabel('radius [mm]'); ylabel('temperature [\circC]');
title('Radial temperature profiles (optimal sleeve position)');
legend(legendEntries, 'Location', 'northwest', 'FontSize', 7);
grid on;
print(f, fullfile(outDir, 'fig3_temperature_profiles.png'), '-dpng', '-r150');
close(f);

% ---- fig4: sleeve temperature + pulsed power (baseline, centered) -------
f = figure('Visible', 'off');
subplot(2,1,1);
plot(Rbase.t/60, Rbase.Tsleeve, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.4); hold on;
plot(xlim(), [P.Tm P.Tm], 'r--', 'LineWidth', 1);
ylabel('sleeve avg. temp [\circC]');
title('Sleeve pulse-heating behaviour (centered case)');
grid on;
subplot(2,1,2);
stairs(Rbase.t/60, Rbase.power, 'Color', [0.27 0.51 0.71], 'LineWidth', 1.2);
xlabel('time [min]'); ylabel('sleeve power [W/m]');
grid on;
print(f, fullfile(outDir, 'fig4_sleeve_pulsing.png'), '-dpng', '-r150');
close(f);

end

% ------------------------------------------------------------------------
function yline_manual(y)
xl = xlim();
plot(xl, [y y], ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.7);
end

function xline_manual(x, color)
yl = ylim();
plot([x x], yl, ':', 'Color', color, 'LineWidth', 1.2);
end


% ------------------------------------------------------------------------
% source: makePlotsFastMelt.m
% ------------------------------------------------------------------------
function makePlotsFastMelt(P, R, Rbare, rowsNS, rowsS, chosen, band, outDir)
%MAKEPLOTSFASTMELT  Result plots for mainFastMelt.m -- mirrors
%PCM_python/fast_melt_design.py's make_plots().

% ---- fig1: gap-width sweeps (no-sleeve time, with-sleeve power) --------
f = figure('Visible', 'off', 'Position', [100 100 1000 420]);

subplot(1,2,1);
gapsNS = [rowsNS.gapMm];
tNS = [rowsNS.t999];
plot(gapsNS, tNS, 'o-', 'Color', [0.27 0.51 0.71], 'LineWidth', 1.4); hold on;
plot(xlim(), [180 180], 'r--', 'LineWidth', 1);
xlabel('PCM annulus gap width [mm]'); ylabel('full-melt time [s]');
title('Pure pipe conduction (no sleeve)');
legend('t999', '3 min target', 'Location', 'northwest');
grid on;

subplot(1,2,2);
gapsS = [rowsS.gapMm];
pS = [rowsS.Pmin];
plot(gapsS, pS, 'o-', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.4); hold on;
plot([chosen.gapMm chosen.gapMm], ylim(), 'g:', 'LineWidth', 1.5);
xlabel('PCM annulus gap width [mm]');
ylabel('min. sleeve power to melt in <170 s [W/m]');
title('With pulsed copper sleeve');
legend('P_{min}', 'recommended', 'Location', 'northwest');
grid on;

print(f, fullfile(outDir, 'fig1_gap_sweeps.png'), '-dpng', '-r150');
close(f);

% ---- fig2: liquid fraction, recommended design vs. no-sleeve at same gap --
f = figure('Visible', 'off');
plot(R.t, R.favg*100, 'LineWidth', 1.8); hold on;
plot(Rbare.t, Rbare.favg*100, '--', 'LineWidth', 1.8);
plot([180 180], ylim(), 'r:', 'LineWidth', 1);
xlabel('time [s]'); ylabel('PCM average liquid fraction [%]');
title(sprintf('Melting progress at the recommended %d mm gap', chosen.gapMm));
legend('with pulsed sleeve', 'no sleeve (same gap)', '3 min target', 'Location', 'southeast');
grid on;
print(f, fullfile(outDir, 'fig2_liquid_fraction.png'), '-dpng', '-r150');
close(f);

% ---- fig3: temperature profile snapshots ---------------------------------
G = createGeometryRadial(P, chosen.frac);
snaps = R.snapshots;
nSnap = numel(snaps);
pick = round(linspace(1, nSnap, min(7, nSnap)));
f = figure('Visible', 'off'); hold on;
cmap = jet(numel(pick));
legendEntries = {};
for i = 1:numel(pick)
    s = snaps(pick(i));
    plot(G.rc*1000, s.T, 'Color', cmap(i,:), 'LineWidth', 1.4);
    legendEntries{end+1} = sprintf('t=%.0f s', s.t); %#ok<AGROW>
end
plot(xlim(), [P.Tm P.Tm], 'r--', 'LineWidth', 1);
legendEntries{end+1} = 'PCM melt point';
yl = ylim();
patch([band(1) band(2) band(2) band(1)]*1000, [yl(1) yl(1) yl(2) yl(2)], ...
      [1 0.65 0], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
legendEntries{end+1} = 'copper sleeve';
xlabel('radius [mm]'); ylabel('temperature [\circC]');
title('Radial temperature profiles (recommended fast-melt design)');
legend(legendEntries, 'Location', 'northeast', 'FontSize', 7);
grid on;
print(f, fullfile(outDir, 'fig3_temperature_profiles.png'), '-dpng', '-r150');
close(f);

% ---- fig4: sleeve power/temperature trace --------------------------------
f = figure('Visible', 'off');
subplot(2,1,1);
plot(R.t, R.Tsleeve, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.4); hold on;
plot(xlim(), [P.Tm P.Tm], 'r--', 'LineWidth', 0.8);
ylabel('sleeve avg. temp [\circC]');
title('Recommended design: sleeve pulse-heating trace');
legend('sleeve temp', 'PCM melt point', 'Location', 'southeast');
grid on;
subplot(2,1,2);
stairs(R.t, R.power, 'Color', [0.27 0.51 0.71], 'LineWidth', 1.2);
xlabel('time [s]'); ylabel('sleeve power [W/m]');
grid on;
print(f, fullfile(outDir, 'fig4_sleeve_pulsing.png'), '-dpng', '-r150');
close(f);

end


% ==========================================================================
% DRIVER -- runs the complete simulation: the geometry as drawn (Part 1)
% and the re-dimensioned <3-minute redesign (Part 2). Everything above
% this point is function definitions; this is the only part that
% actually executes top to bottom.
% ==========================================================================

outDirBase = fullfile(pwd, 'results_master');
outDirAsDrawn = fullfile(outDirBase, 'asDrawn');
outDirFastMelt = fullfile(outDirBase, 'fastMelt');
if ~exist(outDirAsDrawn, 'dir'); mkdir(outDirAsDrawn); end
if ~exist(outDirFastMelt, 'dir'); mkdir(outDirFastMelt); end

fprintf('##########################################################################\n');
fprintf('# PART 1 of 2 -- the geometry as drawn (32.3 mm pipe->wall gap)\n');
fprintf('##########################################################################\n\n');

outDir = outDirAsDrawn;

P = parametersRadial();
matp = materialStruct(P);

fprintf('=== PCM radial model (MATLAB/Octave): hot pipe + pulsed copper sleeve ===\n');
fprintf('rPipe=%.2f mm, rWall=%.2f mm (gap %.2f mm) [container size ASSUMED]\n', ...
    P.rPipe*1000, P.rWall*1000, (P.rWall-P.rPipe)*1000);
fprintf('Thotwater=%.0f C, Tm=%.0f C, Tinit=%.0f C\n', P.Thotwater, P.Tm, P.Tinit);
fprintf('PCM: rho=%.0f, cp=%.0f J/kgK, ks=%.2f, kl=%.2f W/mK, L=%.0f J/kg [ASSUMED dataset]\n\n', ...
    P.rhoPCM, P.cpPCM, P.kPCMs, P.kPCMl, P.LPCM);

%% ---- 1) Baseline: sleeve centered in the gap, as drawn ------------------
fprintf('--- Baseline: sleeve centered in the gap (as drawn) ---\n');
Gbase = createGeometryRadial(P, 0.5);
tic;
[PminBase, TreachedBase] = findMinPowerRadial(P, Gbase, matp);
Rbase = runFullMeltRadial(P, Gbase, matp, PminBase, P.Tm, 20, true);
fprintf('  sleeve band: %.2f-%.2f mm\n', Gbase.rSleeveIn*1000, Gbase.rSleeveOut*1000);
fprintf('  minimum sleeve peak power Pmin = %.1f W/m (reaches %.1f C by tRamp=%.0f s)\n', ...
    PminBase, TreachedBase, P.tRampTarget);
fprintf('  t50=%s  t90=%s  t99=%s  t999(full melt)=%s\n', ...
    formatTimeOrNARadial(Rbase.t50), formatTimeOrNARadial(Rbase.t90), ...
    formatTimeOrNARadial(Rbase.t99), formatTimeOrNARadial(Rbase.t999));
fprintf('  sleeve heater duty cycle: %.1f%%\n', Rbase.dutyCycle*100);
fprintf('  elapsed wall time: %.1f s\n\n', toc);

%% ---- 2) No-sleeve reference ----------------------------------------------
fprintf('--- Reference: no sleeve, pipe conduction only ---\n');
Gref = createGeometryNoSleeve(P);
Rref = runFullMeltRadial(P, Gref, matp, 0, P.Tm, 20, false);
fprintf('  t50=%s  t90=%s  t99=%s  t999(full melt)=%s\n\n', ...
    formatTimeOrNARadial(Rref.t50), formatTimeOrNARadial(Rref.t90), ...
    formatTimeOrNARadial(Rref.t99), formatTimeOrNARadial(Rref.t999));

%% ---- 3) Sleeve-position sweep --------------------------------------------
fprintf('--- Sleeve-position sweep (coarse) ---\n');
fracsCoarse = linspace(0.15, 0.90, 16);
rowsCoarse = runPositionSweepRadial(P, matp, fracsCoarse, 40, true);

[~, iBest] = min([rowsCoarse.t999]);
bestCoarse = rowsCoarse(iBest);
fprintf('\n  coarse best: frac=%.2f -> t999=%.0f s\n', bestCoarse.frac, bestCoarse.t999);

fprintf('\n--- Sleeve-position sweep (refine near coarse optimum) ---\n');
lo = max(0.10, bestCoarse.frac - 0.08);
hi = min(0.95, bestCoarse.frac + 0.08);
fracsFine = linspace(lo, hi, 9);
rowsFine = runPositionSweepRadial(P, matp, fracsFine, 40, true);

rowsAll = [rowsCoarse, rowsFine];
[~, iBestAll] = min([rowsAll.t999]);
best = rowsAll(iBestAll);

fprintf('\n  OPTIMAL sleeve position: frac=%.3f (rMid=%.2f mm)\n', best.frac, best.rMid*1000);
fprintf('  Pmin at optimum = %.1f W/m\n', best.Pmin);
fprintf('  t99=%s   t999(full melt)=%s\n', formatTimeOrNARadial(best.t99), formatTimeOrNARadial(best.t999));

%% ---- Full run at the optimal position (for plotting) ----------------------
Gopt = createGeometryRadial(P, best.frac);
Ropt = runFullMeltRadial(P, Gopt, matp, best.Pmin, P.Tm, 20, true);

%% ---- Save CSV of the sweep -------------------------------------------------
[~, order] = sort([rowsAll.frac]);
rowsSorted = rowsAll(order);
fid = fopen(fullfile(outDir, 'sweep.csv'), 'w');
fprintf(fid, 'frac,rIn,rOut,rMid,Pmin_W_per_m,t50,t90,t99,t999,dutyCycle,converged\n');
for i = 1:numel(rowsSorted)
    r = rowsSorted(i);
    fprintf(fid, '%.4f,%.6f,%.6f,%.6f,%.4f,%s,%s,%s,%s,%.4f,%d\n', ...
        r.frac, r.rIn, r.rOut, r.rMid, r.Pmin, ...
        numOrEmpty(r.t50), numOrEmpty(r.t90), numOrEmpty(r.t99), numOrEmpty(r.t999), ...
        r.dutyCycle, r.converged);
end
fclose(fid);

%% ---- Plots -----------------------------------------------------------------
makePlotsAsDrawn(P, Rbase, Rref, Ropt, rowsSorted, best, ...
    [Gbase.rSleeveIn, Gbase.rSleeveOut], [Gopt.rSleeveIn, Gopt.rSleeveOut], outDir);

%% ---- Summary -----------------------------------------------------------------
fid = fopen(fullfile(outDir, 'summary.txt'), 'w');
fprintf(fid, 'rPipe_mm=%.4f\nrWall_mm=%.4f\nThotwater=%.2f\nTm=%.2f\nTinit=%.2f\ntRampTarget=%.1f\n\n', ...
    P.rPipe*1000, P.rWall*1000, P.Thotwater, P.Tm, P.Tinit, P.tRampTarget);
fprintf(fid, '[baseline_centered]\nPmin_W_per_m=%.4f\nt50=%s\nt90=%s\nt99=%s\nt999=%s\ndutyCycle=%.4f\n\n', ...
    PminBase, numOrEmpty(Rbase.t50), numOrEmpty(Rbase.t90), numOrEmpty(Rbase.t99), numOrEmpty(Rbase.t999), Rbase.dutyCycle);
fprintf(fid, '[no_sleeve_reference]\nt50=%s\nt90=%s\nt99=%s\nt999=%s\n\n', ...
    numOrEmpty(Rref.t50), numOrEmpty(Rref.t90), numOrEmpty(Rref.t99), numOrEmpty(Rref.t999));
fprintf(fid, '[optimal]\nfrac=%.4f\nrMid_mm=%.4f\nPmin_W_per_m=%.4f\nt99=%s\nt999=%s\n', ...
    best.frac, best.rMid*1000, best.Pmin, numOrEmpty(best.t99), numOrEmpty(best.t999));
fclose(fid);

fprintf('\nWrote plots, sweep.csv and summary.txt to %s/\n', outDir);
fprintf('\n=== PART 1 HEADLINE RESULT ===\n');
fprintf('As drawn (sleeve centered, %.0f W/m minimum pulse power): full melt in %s\n', ...
    PminBase, formatTimeOrNARadial(Rbase.t999));
fprintf('No auxiliary sleeve heater at all (pipe conduction only): full melt in %s\n', ...
    formatTimeOrNARadial(Rref.t999));
fprintf('Minimum achievable melting time (sleeve at %.0f%% pipe->wall, %.0f W/m): full melt in %s\n\n', ...
    best.frac*100, best.Pmin, formatTimeOrNARadial(best.t999));


fprintf('##########################################################################\n');
fprintf('# PART 2 of 2 -- re-dimensioned for melting in under 3 minutes\n');
fprintf('##########################################################################\n\n');

outDir = outDirFastMelt;

TARGET_S = 180.0;         % the <3-minute requirement
BISECT_TARGET_S = 170.0;  % bisect inside the target for margin
T_SAFETY_CAP = 180.0;     % deg C, sleeve heater cuts off above this
REALISTIC_POWER_CAP = 2000.0;  % W/m -- "realistic" cutoff for the recommended design

Pbase = parametersRadial();
rPipe = Pbase.rPipe;

fprintf('=== Fast-melt design: re-dimension for <3-minute melting ===\n\n');

%% ---- Step 1: pure pipe conduction (no sleeve) --------------------------
fprintf('--- Step 1: pure pipe conduction (no sleeve) -- gap-width sweep ---\n');
gapsNoSleeve = [2 3 4 5 6 8 10 12];
rowsNS = sweepNoSleeveGaps(rPipe, gapsNoSleeve, 220, 0.25, 1200, Pbase, true);

%% ---- Step 2: with pulsed copper sleeve ----------------------------------
fprintf('\n--- Step 2: with pulsed copper sleeve -- gap-width sweep (target <=%.0f s, cap %.0f C) ---\n', ...
    BISECT_TARGET_S, T_SAFETY_CAP);
gapsSleeve = [6 7 8 9 10 11 12 13 14];
rowsS = sweepWithSleeveGaps(rPipe, gapsSleeve, 260, 0.2, 600, 0.001, Pbase, ...
    BISECT_TARGET_S, T_SAFETY_CAP, true);

%% ---- Pick the recommended design ----------------------------------------
feasible = rowsS(~isnan([rowsS.Pmin]));
candidates = feasible([feasible.Pmin] <= REALISTIC_POWER_CAP);
[~, iMax] = max([candidates.gapMm]);
chosen = candidates(iMax);
fprintf('\nRecommended design: gap=%d mm, frac=%.2f, Pmin=%.1f W/m -> t999=%.1f s\n', ...
    chosen.gapMm, chosen.frac, chosen.Pmin, chosen.t999);

gapMm = chosen.gapMm;
frac = chosen.frac;
Pdesign = round((chosen.Pmin * 1.05)/10)*10;  % +5% margin, rounded to 10 W/m

P = Pbase;
P.rPipe = rPipe; P.rWall = rPipe + gapMm/1000;
P.nCells = 400; P.dt = 0.1; P.tMax = 400; P.tSleeve = 0.001;
matp = materialStruct(P);

G = createGeometryRadial(P, frac);
R = runFullMeltRadial(P, G, matp, Pdesign, T_SAFETY_CAP, 20, true);
fprintf('Final verified run (finer grid, +5%% power margin = %.0f W/m): t50=%.1fs t90=%.1fs t99=%.1fs t999=%.1fs\n', ...
    Pdesign, R.t50, R.t90, R.t99, R.t999);

%% ---- No-sleeve check at the same gap, for contrast ------------------------
Gbare = createGeometryNoSleeve(P);
Rbare = runFullMeltRadial(P, Gbare, matp, 0, P.Tm, 20, false);

%% ---- Save CSVs --------------------------------------------------------------
fid = fopen(fullfile(outDir, 'no_sleeve_sweep.csv'), 'w');
fprintf(fid, 'gap_mm,t999\n');
for i = 1:numel(rowsNS)
    fprintf(fid, '%d,%s\n', rowsNS(i).gapMm, numOrEmpty(rowsNS(i).t999));
end
fclose(fid);

fid = fopen(fullfile(outDir, 'sleeve_sweep.csv'), 'w');
fprintf(fid, 'gap_mm,frac,P_min_W_per_m,t999\n');
for i = 1:numel(rowsS)
    fprintf(fid, '%d,%.4f,%s,%s\n', rowsS(i).gapMm, rowsS(i).frac, ...
        numOrEmpty(rowsS(i).Pmin), numOrEmpty(rowsS(i).t999));
end
fclose(fid);

fid = fopen(fullfile(outDir, 'summary.txt'), 'w');
fprintf(fid, 'target_s=%.1f\n\n[recommended]\nr_pipe_mm=%.4f\ngap_mm=%d\nr_wall_mm=%.4f\n', ...
    TARGET_S, P.rPipe*1000, gapMm, P.rWall*1000);
fprintf(fid, 'sleeve_thickness_mm=%.4f\nsleeve_center_frac=%.4f\nsleeve_band_mm=%.4f-%.4f\n', ...
    P.tSleeve*1000, frac, G.rSleeveIn*1000, G.rSleeveOut*1000);
fprintf(fid, 'P_design_W_per_m=%.1f\nT_safety_cap_C=%.1f\n', Pdesign, T_SAFETY_CAP);
fprintf(fid, 't50=%.2f\nt90=%.2f\nt99=%.2f\nt999=%.2f\n', R.t50, R.t90, R.t99, R.t999);
fclose(fid);

%% ---- Plots --------------------------------------------------------------
makePlotsFastMelt(P, R, Rbare, rowsNS, rowsS, chosen, [G.rSleeveIn, G.rSleeveOut], outDir);

fprintf('\nWrote plots, CSVs and summary.txt to %s/\n', outDir);
fprintf('\n=== PART 2 FINAL RECOMMENDATION ===\n');
fprintf('Pipe radius (fixed, given): %.2f mm\n', P.rPipe*1000);
fprintf('Container inner radius: %.2f mm (PCM annulus thickness: %d mm, was 32.3 mm)\n', ...
    P.rWall*1000, gapMm);
fprintf('Copper sleeve: %.1f mm thick, centered at %.0f%% of the way from pipe to wall (%.2f-%.2f mm)\n', ...
    P.tSleeve*1000, frac*100, G.rSleeveIn*1000, G.rSleeveOut*1000);
fprintf('Sleeve pulse power: %.0f W/m of pipe length (e.g. %.0f W for a 0.3 m section), capped off above %.0f C\n', ...
    Pdesign, Pdesign*0.3, T_SAFETY_CAP);
fprintf('Result: fully melted (99.9%%) in %.0f s = %.2f min (target: < 180 s)\n', R.t999, R.t999/60);
fprintf('Without the sleeve at this same gap (400 s cutoff): t50=%s t90=%s t99=%s t999=%s\n\n', ...
    numOrEmpty(Rbare.t50), numOrEmpty(Rbare.t90), numOrEmpty(Rbare.t99), numOrEmpty(Rbare.t999));

fprintf('##########################################################################\n');
fprintf('# COMPLETE. Results written under: %s\n', outDirBase);
fprintf('#   %s  (Part 1: geometry as drawn)\n', outDirAsDrawn);
fprintf('#   %s  (Part 2: <3-minute redesign)\n', outDirFastMelt);
fprintf('##########################################################################\n');
