function [fl, dfldT] = liquidFraction(T, Ts, Tl)
%LIQUIDFRACTION  Smooth (cubic Hermite / "smoothstep") liquid fraction
%and its temperature derivative, used for the apparent-heat-capacity
%phase-change model.
%
%   [fl, dfldT] = LIQUIDFRACTION(T, Ts, Tl)
%
%       fl(T) = 0                          , T <= Ts
%              = 3*xi^2 - 2*xi^3            , Ts < T < Tl
%              = 1                          , T >= Tl
%   where xi = (T - Ts) / (Tl - Ts).
%
%   dfldT(T) = 6*xi*(1-xi) / (Tl - Ts)  inside the mushy interval and 0
%   outside it.
%
%   T may be any array; fl and dfldT are returned the same size/shape.

dT = Tl - Ts;
if dT <= 0
    error('liquidFraction:badInterval', 'Require Tl > Ts.');
end

xi = (T - Ts) / dT;
xi = min(max(xi, 0), 1);          % clamp to [0,1]

fl    = 3*xi.^2 - 2*xi.^3;
dfldT = 6*xi.*(1 - xi) / dT;

% Outside the mushy interval fl is exactly flat (0 or 1), so the
% derivative must be exactly zero there even though the clamped xi
% formula above already gives 0 at xi=0 and xi=1 -- kept explicit for
% clarity/robustness against round-off.
outside = (T <= Ts) | (T >= Tl);
dfldT(outside) = 0;

end
