function [t999, R] = t999ForPower(P, G, matp, power, Tsetpoint, recordEvery)
%T999FORPOWER  Convenience wrapper: run runFullMeltRadial at a given
%sleeve peak power (bang-bang up to Tsetpoint) and return just the
%full-melt (99.9%) milestone, plus the full result struct.

if nargin < 6 || isempty(recordEvery), recordEvery = 100; end
R = runFullMeltRadial(P, G, matp, power, Tsetpoint, recordEvery, false);
t999 = R.t999;

end
