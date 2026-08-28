function [Tnew, fl, k, Capp, nPicard] = solveOneTimeStep(Tn, Qv, G, C, P)
%SOLVEONETIMESTEP  Advance the temperature field by one implicit
%(backward-Euler) time step, solving
%
%   Capp(T) * V/dt * (T - Tn) = sum_faces G_face*(T_nb - T) + Q'''*V
%
%on the reduced active-cell grid G. Because k(T) and Capp(T) are
%nonlinear (through the liquid fraction), a few Picard sub-iterations
%are used: the matrix is assembled with properties evaluated at the
%current temperature guess, solved, and the guess updated, until the
%change in T is below P.picardTol or P.picardIters is reached.
%
%   [Tnew, fl, k, Capp, nPicard] = SOLVEONETIMESTEP(Tn, Qv, G, C, P)
%
%       Tn   : Na x 1 temperature at the start of the step [deg C]
%       Qv   : Na x 1 volumetric heat generation            [W/m^3]
%       G    : geometry struct from createGeometry.m
%       C    : composite-property struct from compositeProperties.m
%       P    : parameter struct from parameters.m
%
%       Tnew : Na x 1 temperature at the end of the step
%       fl   : Na x 1 liquid fraction (0 for heater cells)
%       k    : Na x 1 conductivity used in the final iteration
%       Capp : Na x 1 apparent heat capacity used in the final iteration
%       nPicard : number of Picard iterations actually performed

Na  = G.Na;
pcm = G.pcmMaskA;
heater = ~pcm;               % heater1 + heater2 cells
dt  = P.dt;
V   = G.V;

Tguess = Tn;

for it = 1:P.picardIters
    fl = liquidFraction(Tguess, P.Ts, P.Tl);
    fl(heater) = 0;          % heater material never undergoes PCM phase change

    k = zeros(Na, 1);
    k(pcm)    = (1 - fl(pcm)).*C.ks_eff + fl(pcm).*C.kl_eff;
    k(heater) = C.k_Cu;

    % Energy-conserving ("secant") apparent heat capacity -- see
    % apparentCapacity.m for why the naive tangent-derivative form is
    % not used here.
    Capp = zeros(Na, 1);
    Capp(pcm)    = apparentCapacity(Tguess(pcm), Tn(pcm), P.Ts, P.Tl, ...
                                     C.rhocp_eff, C.rhoL_eff);
    Capp(heater) = C.rhocp_Cu;

    A = assembleConductionMatrix(G, k, Capp, V, dt);
    b = Capp .* (V/dt) .* Tn + Qv * V;

    Tnew = A \ b;

    change = norm(Tnew - Tguess, Inf);
    scale  = max(1, norm(Tguess, Inf));
    Tguess = Tnew;

    if change < P.picardTol * scale
        break;
    end
end

nPicard = it;

% Report fl/k/Capp consistent with the returned temperature.
[fl, ~] = liquidFraction(Tnew, P.Ts, P.Tl);
fl(heater) = 0;

end

% ------------------------------------------------------------------------
function A = assembleConductionMatrix(G, k, Capp, V, dt)
%ASSEMBLECONDUCTIONMATRIX  Sparse FV matrix for cubic cells of size h.
%Face conductance between adjacent active cells p,q is the harmonic
%mean of k(p),k(q) times h (= face area h^2 / centre spacing h). Cells
%with no active neighbour on a face are implicitly adiabatic there
%(no entry is added for that face), which reproduces the insulated
%outer-cylinder / box boundary condition.

h  = G.h;
Na = G.Na;

pairsAll = [G.pairs.x; G.pairs.y; G.pairs.z];
p = pairsAll(:,1);
q = pairsAll(:,2);

kH = 2 * k(p) .* k(q) ./ (k(p) + k(q));   % harmonic mean conductivity
Gc = kH * h;                              % face conductance [W/K]

rows = [p; q; p; q];
cols = [q; p; p; q];
vals = [-Gc; -Gc; Gc; Gc];

% Time-derivative (capacity) diagonal
rows = [rows; (1:Na)'];
cols = [cols; (1:Na)'];
vals = [vals; Capp .* (V/dt)];

A = sparse(rows, cols, vals, Na, Na);

end
