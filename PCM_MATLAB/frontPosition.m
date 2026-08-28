function xf = frontPosition(xc, T, Tm, tol)
%FRONTPOSITION  Linear-interpolate the x location where a 1-D
%temperature profile T(xc) crosses the melting temperature Tm
%(T assumed decreasing with x, hot wall at x=0). Used by
%validateStefan.m to extract the numerical melt-front position.
%
%   In the one-phase Stefan setup the *unheated* far-field is exactly
%   at T = Tm (not below it), so a non-strict T>=Tm test would flag
%   the whole (still solid) domain as melted from t=0. TOL (default
%   1e-6 deg C) enforces a strict, numerically-safe "has actually
%   warmed above Tm" test.
%
%   Returns xc(1) if the whole domain is still solid (should not
%   normally happen since the wall cell is always warmest) and xc(end)
%   if the whole domain has melted.

if nargin < 4
    tol = 1e-6;
end

idx = find(T > Tm + tol, 1, 'last');
if isempty(idx)
    xf = xc(1);
elseif idx == numel(T)
    xf = xc(end);
else
    T1 = T(idx); T2 = T(idx+1);
    x1 = xc(idx); x2 = xc(idx+1);
    frac = (Tm - T1) / (T2 - T1);
    xf = x1 + frac*(x2 - x1);
end

end
