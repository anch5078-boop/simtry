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
