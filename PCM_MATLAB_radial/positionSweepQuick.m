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
