function [fl, dfldT] = liquidFraction(T, Ts, Tl)
%LIQUIDFRACTION  Smooth (cubic Hermite / "smoothstep") liquid fraction
%and its temperature derivative -- identical formula to
%PCM_MATLAB/liquidFraction.m, used here by Model B's apparent-heat-
%capacity phase-change model.
%
%   [fl, dfldT] = LIQUIDFRACTION(T, Ts, Tl)
%
%       fl(T) = 0                          , T <= Ts
%              = 3*xi^2 - 2*xi^3            , Ts < T < Tl
%              = 1                          , T >= Tl
%   where xi = (T - Ts) / (Tl - Ts).
%
%   T may be any array; fl and dfldT are returned the same size/shape.

dT = Tl - Ts;
if dT <= 0
    error('liquidFraction:badInterval', 'Require Tl > Ts.');
end

xi = (T - Ts) / dT;
xi = min(max(xi, 0), 1);

fl    = 3*xi.^2 - 2*xi.^3;
dfldT = 6*xi.*(1 - xi) / dT;

outside = (T <= Ts) | (T >= Tl);
dfldT(outside) = 0;

end
