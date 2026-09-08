function d3 = deltaEffCubed(delta, b)
%DELTAEFFCUBED  delta_eff^3 = delta^3*(delta+4b)/(delta+b) -- the
%slip-modified effective hydrodynamic gap (cubed) for one-wall
%Navier-slip plane Poiseuille flow, derived in ccmSlipModel.m's header
%comment. Reduces to delta.^3 at b=0.
%
%   d3 = DELTAEFFCUBED(delta, b)

d3 = delta.^3 .* (delta + 4*b) ./ (delta + b);

end
