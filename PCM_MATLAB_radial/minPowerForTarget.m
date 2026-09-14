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
