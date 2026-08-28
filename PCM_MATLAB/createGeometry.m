function G = createGeometry(P)
%CREATEGEOMETRY  Build the structured Cartesian grid, cylindrical PCM
%mask and the two cylindrical-shell heater masks.
%
%   G = CREATEGEOMETRY(P) returns a struct describing the finite-volume
%   grid: only cells that satisfy x^2+y^2 <= R^2 (and 0<=z<=H) are
%   "active" -- everything else is simply not simulated, which is
%   equivalent to a perfectly insulated (adiabatic) outer cylindrical
%   wall. Active cells are further split into PCM cells and heater-1 /
%   heater-2 cells.
%
%   To keep the time-step solver fast, only ACTIVE cells are carried as
%   unknowns (T is a column vector of length G.Na, not a full 3-D
%   array). G.activeIdx / G.g2a map between the full nx*ny*nz grid and
%   this reduced active-cell numbering, and G.pairs.{x,y,z} list the
%   face-adjacent pairs of active cells (already reduced-numbered) so
%   solveOneTimeStep.m can assemble the sparse conduction matrix
%   without ever touching inactive cells.

h = P.h;

% Cell-centred coordinate vectors. The domain box exactly bounds the
% cylinder: x,y in [-R,R], z in [0,H].
nx = round(2*P.R/h);
ny = nx;
nz = round(P.H/h);

xc = linspace(-P.R + h/2, P.R - h/2, nx);
yc = linspace(-P.R + h/2, P.R - h/2, ny);
zc = linspace(h/2,        P.H - h/2, nz);

[X, Y, Z] = ndgrid(xc, yc, zc);
r = sqrt(X.^2 + Y.^2);

%% ---- Masks (full nx x ny x nz logical arrays) ------------------------
cylinderMask = r <= P.R;

heater1Mask = cylinderMask & r >= P.r1i & r <= P.r1o & ...
              Z >= P.zbase & Z <= P.zbase + P.Hheater;

heater2Mask = cylinderMask & r >= P.r2i & r <= P.r2o & ...
              Z >= P.zbase & Z <= P.zbase + P.Hheater;

pcmMask   = cylinderMask & ~heater1Mask & ~heater2Mask;
activeMask = cylinderMask;   % pcm + heater1 + heater2

if ~any(pcmMask(:))
    error('createGeometry:noPCM', 'No PCM cells found -- check geometry.');
end
if ~any(heater1Mask(:)) || ~any(heater2Mask(:))
    error('createGeometry:noHeater', ...
        'Heater mask is empty -- check radii against grid resolution h.');
end

%% ---- Reduced active-cell numbering ------------------------------------
N  = nx*ny*nz;
L  = reshape(1:N, nx, ny, nz);      % global linear index of every cell

activeIdx = find(activeMask(:));
Na = numel(activeIdx);

g2a = zeros(N, 1);                  % global idx -> active idx (0 = inactive)
g2a(activeIdx) = (1:Na)';

matType = zeros(Na, 1);             % 0 = PCM, 1 = heater1, 2 = heater2
matType(g2a(heater1Mask(:))) = 1;
matType(g2a(heater2Mask(:))) = 2;

%% ---- Face-adjacent active-cell pairs (x, y, z directions) -------------
G.pairs.x = neighborPairs(activeMask, L, g2a, 1);
G.pairs.y = neighborPairs(activeMask, L, g2a, 2);
G.pairs.z = neighborPairs(activeMask, L, g2a, 3);

%% ---- Volumes -----------------------------------------------------------
V   = h^3;                                  % volume of one cubic cell
Vh1 = nnz(heater1Mask) * h^3;                % heater-1 total volume
Vh2 = nnz(heater2Mask) * h^3;                % heater-2 total volume

%% ---- Pack output ---------------------------------------------------------
G.h  = h;
G.nx = nx; G.ny = ny; G.nz = nz; G.N = N;
G.xc = xc; G.yc = yc; G.zc = zc;
G.X = X; G.Y = Y; G.Z = Z; G.r = r;

G.cylinderMask = cylinderMask;
G.pcmMask      = pcmMask;
G.heater1Mask  = heater1Mask;
G.heater2Mask  = heater2Mask;
G.activeMask   = activeMask;

G.activeIdx = activeIdx;
G.g2a       = g2a;
G.matType   = matType;              % Na x 1, per ACTIVE cell
G.pcmMaskA     = matType == 0;      % logical, Na x 1
G.heater1MaskA = matType == 1;
G.heater2MaskA = matType == 2;

G.Na = Na;
G.V   = V;
G.Vh1 = Vh1;
G.Vh2 = Vh2;

end

% ------------------------------------------------------------------------
function pairs = neighborPairs(activeMask, L, g2a, dim)
%NEIGHBORPAIRS  Reduced-numbered [p q] pairs of face-adjacent active
%cells along dimension DIM (1=x, 2=y, 3=z).

switch dim
    case 1
        maskA = activeMask(1:end-1, :, :);
        maskB = activeMask(2:end,   :, :);
        idxA  = L(1:end-1, :, :);
        idxB  = L(2:end,   :, :);
    case 2
        maskA = activeMask(:, 1:end-1, :);
        maskB = activeMask(:, 2:end,   :);
        idxA  = L(:, 1:end-1, :);
        idxB  = L(:, 2:end,   :);
    case 3
        maskA = activeMask(:, :, 1:end-1);
        maskB = activeMask(:, :, 2:end);
        idxA  = L(:, :, 1:end-1);
        idxB  = L(:, :, 2:end);
end

pairMask = maskA & maskB;
p = g2a(idxA(pairMask));
q = g2a(idxB(pairMask));
pairs = [p(:) q(:)];

end
