function [rIn, rOut] = sleeveRadii(P, centerFrac)
%SLEEVERADII  Inner/outer radius of the copper sleeve for a given
%fractional position across the pipe->wall gap (0 = touching the pipe,
%1 = touching the wall). Defaults to P.centerFrac ("centered in the
%gap", as drawn) when centerFrac is omitted.

if nargin < 2 || isempty(centerFrac)
    centerFrac = P.centerFrac;
end

rMid = P.rPipe + centerFrac * (P.rWall - P.rPipe);
rIn  = rMid - P.tSleeve/2;
rOut = rMid + P.tSleeve/2;

end
