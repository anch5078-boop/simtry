function Ffull = expandToFullGrid(fieldActive, G, fillValue)
%EXPANDTOFULLGRID  Scatter a reduced active-cell field (Na x 1) back
%onto the full nx-by-ny-by-nz Cartesian grid for slicing/plotting.
%Inactive cells (outside the PCM cylinder) are set to FILLVALUE
%(default NaN, so they simply don't plot).
%
%   Ffull = EXPANDTOFULLGRID(fieldActive, G, fillValue)

if nargin < 3
    fillValue = NaN;
end

Ffull = fillValue * ones(G.nx, G.ny, G.nz);
Ffull(G.activeIdx) = fieldActive;

end
