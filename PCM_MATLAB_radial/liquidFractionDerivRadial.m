function dfldT = liquidFractionDerivRadial(T, Ts, Tl)
%LIQUIDFRACTIONDERIVRADIAL  d(fl)/dT for the smoothed liquid fraction in
%liquidFractionRadial.m -- used only as the tiny-step fallback inside
%apparentCapacityRadial.m.

dT = Tl - Ts;
xi = min(max((T - Ts) / dT, 0), 1);
dfldT = 6*xi.*(1 - xi) / dT;

outside = (T <= Ts) | (T >= Tl);
dfldT(outside) = 0;

end
