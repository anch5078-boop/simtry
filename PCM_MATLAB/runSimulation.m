function R = runSimulation(P, G, C, verbose)
%RUNSIMULATION  Main Version-1 time-marching loop:
%
%   Initialize T
%     -> identify PCM / heater cells            (done in createGeometry)
%     -> at each time step:
%          check pulse ON/OFF & apply Q''' in heater cells
%          compute liquid fraction, k(T), Capp(T)
%          solve implicit heat equation
%          update T
%          compute average PCM liquid fraction
%          check t50 / t90 / t99
%     -> repeat until t99 reached or t reaches P.tMax
%
%   R = RUNSIMULATION(P, G, C, verbose)
%
%       P, G, C  : from parameters.m, createGeometry.m, compositeProperties.m
%       verbose  : (optional, default true) print progress to console
%
%   Returns a results struct R with time histories, t50/t90/t99, and
%   periodic full-field snapshots for plotting.

if nargin < 4
    verbose = true;
end

Na = G.Na;
T  = P.Tinit * ones(Na, 1);
t  = 0;

nSteps = ceil(P.tMax / P.dt);

tHist          = zeros(nSteps+1, 1);
favgHist       = zeros(nSteps+1, 1);
TheaterMaxHist = zeros(nSteps+1, 1);
TmaxHist       = zeros(nSteps+1, 1);
PtotHist       = zeros(nSteps+1, 1);

t50 = NaN; t90 = NaN; t99 = NaN;

snapshots = struct('t', {}, 'T', {}, 'fl', {});

% ---- record the initial state (step 0) ----
[fl0, ~] = liquidFraction(T, P.Ts, P.Tl);
fl0(~G.pcmMaskA) = 0;
tHist(1)          = 0;
favgHist(1)        = mean(fl0(G.pcmMaskA));
TheaterMaxHist(1)  = max(T(~G.pcmMaskA));
TmaxHist(1)        = max(T);
PtotHist(1)        = 0;
snapshots(end+1) = struct('t', 0, 'T', T, 'fl', fl0);

nRec = 1;
for step = 1:nSteps
    tNew = t + P.dt;

    [Qv, P1, P2] = heaterPower(tNew, G, P); %#ok<ASGLU>

    [Tnew, fl, ~, ~, ~] = solveOneTimeStep(T, Qv, G, C, P);

    favg = mean(fl(G.pcmMaskA));   % equal cell volumes -> plain mean

    nRec = nRec + 1;
    tHist(nRec)          = tNew;
    favgHist(nRec)        = favg;
    TheaterMaxHist(nRec)  = max(Tnew(~G.pcmMaskA));
    TmaxHist(nRec)        = max(Tnew);
    PtotHist(nRec)        = P1 + P2;

    if isnan(t50) && favg >= 0.50
        t50 = tNew;
        if verbose, fprintf('  t50 reached at t = %.2f s\n', t50); end
    end
    if isnan(t90) && favg >= 0.90
        t90 = tNew;
        if verbose, fprintf('  t90 reached at t = %.2f s\n', t90); end
    end
    if isnan(t99) && favg >= 0.99
        t99 = tNew;
        if verbose, fprintf('  t99 reached at t = %.2f s\n', t99); end
    end

    if mod(step, P.saveEvery) == 0 || ~isnan(t99)
        snapshots(end+1) = struct('t', tNew, 'T', Tnew, 'fl', fl); %#ok<AGROW>
    end

    if verbose && mod(step, max(1, round(nSteps/20))) == 0
        fprintf('t = %7.1f s   favg = %5.3f   Tmax = %6.2f C   Theater,max = %6.2f C\n', ...
            tNew, favg, max(Tnew), max(Tnew(~G.pcmMaskA)));
    end

    T = Tnew;
    t = tNew;

    if ~isnan(t99)
        break;   % objective reached -- no need to keep integrating
    end
end

% ---- trim preallocated arrays to what was actually used ----
tHist          = tHist(1:nRec);
favgHist       = favgHist(1:nRec);
TheaterMaxHist = TheaterMaxHist(1:nRec);
TmaxHist       = TmaxHist(1:nRec);
PtotHist       = PtotHist(1:nRec);

if isnan(t99) && verbose
    fprintf(['WARNING: t99 not reached within tMax = %.1f s ', ...
             '(final favg = %.3f). Increase P.tMax.\n'], P.tMax, favgHist(end));
end

R.t          = tHist;
R.favg       = favgHist;
R.TheaterMax = TheaterMaxHist;
R.Tmax       = TmaxHist;
R.Ptot       = PtotHist;
R.t50 = t50; R.t90 = t90; R.t99 = t99;
R.snapshots  = snapshots;
R.Tfinal     = T;
R.P = P;

end
