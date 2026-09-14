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
