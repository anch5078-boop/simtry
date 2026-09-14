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
