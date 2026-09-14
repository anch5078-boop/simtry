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
