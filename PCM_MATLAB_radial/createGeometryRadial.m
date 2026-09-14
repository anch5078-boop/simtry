function G = createGeometryRadial(P, centerFrac)
%CREATEGEOMETRYRADIAL  Build the uniform-spacing 1-D radial finite-volume
%grid from P.rPipe to P.rWall, and the copper-sleeve cell mask.
%
%   G = CREATEGEOMETRYRADIAL(P, centerFrac)
%
%       P          : parameter struct from parametersRadial.m
%       centerFrac : (optional) sleeve mid-radius as a fraction of the
%                    pipe->wall gap; defaults to P.centerFrac.
%
%   All of G's "volume"/"area" fields are per unit AXIAL length of pipe
%   (the standard reduction of an axisymmetric, axially-uniform problem
%   to one radial dimension): V(i) is really an area [m^2], Aface(i) is
%   really a circumference [m].

if nargin < 2
    centerFrac = P.centerFrac;
end

N = P.nCells;
rFace = linspace(P.rPipe, P.rWall, N+1)';
rc    = 0.5*(rFace(1:end-1) + rFace(2:end));
dr    = rFace(2) - rFace(1);
V     = pi*(rFace(2:end).^2 - rFace(1:end-1).^2);
Aface = 2*pi*rFace;

[rSleeveIn, rSleeveOut] = sleeveRadii(P, centerFrac);
isCu = rc >= rSleeveIn & rc <= rSleeveOut;

if ~any(isCu)
    error('createGeometryRadial:noSleeve', ...
        'No grid cell falls inside the sleeve band [%.5g, %.5g] m -- increase nCells or tSleeve.', ...
        rSleeveIn, rSleeveOut);
end

G.N = N;
G.rFace = rFace;
G.rc = rc;
G.dr = dr;
G.V = V;
G.Aface = Aface;
G.isCu = isCu;
G.isPcm = ~isCu;
G.rSleeveIn = rSleeveIn;
G.rSleeveOut = rSleeveOut;

end
