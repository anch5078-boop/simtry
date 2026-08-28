%VALIDATESTEFAN  Validate the apparent-heat-capacity conduction/phase-
%change scheme against the classical 1-D one-phase Stefan melting
%problem, for which an exact analytical melt-front position exists.
%
%   Semi-infinite PCM, uniformly at the melting temperature Tm,
%   suddenly exposed at x=0 to a fixed wall temperature T0 > Tm.
%   Because the solid is initially (and remains) exactly at Tm, it
%   carries no sensible heat and the exact 1-D "Neumann/Stefan"
%   solution for the melt-front position is
%
%       s(t) = 2*lambda*sqrt(alpha*t)
%
%   where alpha = k/(rho*cp) (liquid-phase diffusivity) and lambda
%   solves
%
%       lambda*exp(lambda^2)*erf(lambda) = St / sqrt(pi),
%       St = cp*(T0-Tm) / L        (Stefan number)
%
%   This script builds a plain 1-D implicit finite-volume solver using
%   the same apparent-heat-capacity approach (liquidFraction.m,
%   apparentCapacity.m) as the 3-D solver, runs it, extracts the
%   numerical melt-front position by locating the T = Tm isotherm, and
%   reports the error against s(t).
%
%   NOTE on the initial condition: the domain is initialized at Ts
%   (the bottom of the smoothed mushy interval, f_l=0 exactly), not
%   exactly at Tm. Initializing exactly at Tm would place every
%   untouched "far field" cell precisely at the midpoint of the
%   smoothed liquid-fraction curve (f_l(Tm)=0.5 by construction) --
%   purely an artifact of the smoothing, not real absorbed latent
%   heat -- which biases a "T crosses Tm" front diagnostic by several
%   tens of percent even though the scheme conserves energy exactly.
%   Starting at Ts (indistinguishable from Tm once dTmush is small)
%   removes that ambiguity and lets the T=Tm isotherm crossing serve
%   as a clean measure of the true melt front.

clear; clc; close all;
addpath(pwd);

fprintf('=== 1-D Stefan problem validation ===\n');

%% ---- Problem definition ------------------------------------------------
k    = 0.20;     % W/(m K)
rho  = 800;      % kg/m^3
cp   = 2000;     % J/(kg K)
L    = 180000;   % J/kg
Tm   = 60;       % deg C, melting point
T0   = 80;       % deg C, wall (Dirichlet) temperature

% Smoothed mushy interval used by the numerical (enthalpy) scheme.
% As dTmush -> 0 the numerical scheme -> the sharp-interface problem.
dTmush = 0.4;    % deg C
Ts = Tm - dTmush/2;
Tl = Tm + dTmush/2;

Tinf = Ts;       % initial / far-field temperature: f_l=0 exactly (see
                 % NOTE above) -- physically indistinguishable from the
                 % one-phase assumption's "solid uniformly at Tm" once
                 % dTmush is small.

alpha = k/(rho*cp);
St    = cp*(T0-Tm)/L;

lambda = fzero(@(lam) lam.*exp(lam.^2).*erf(lam) - St/sqrt(pi), 0.5);
fprintf('Stefan number St = %.4f, lambda = %.5f\n', St, lambda);

%% ---- Numerical grid / time -----------------------------------------
tEnd = 600;                              % s
sEnd = 2*lambda*sqrt(alpha*tEnd);
Lx   = 3*sEnd;                           % domain long enough vs. front travel

nx = 500;
dx = Lx/nx;
xc = dx/2 + (0:nx-1)'*dx;                % cell centers

dt = 0.05;
nSteps = round(tEnd/dt);

Gc_int = k/dx;        % interior face conductance (per unit area), constant k
Gc_bc  = 2*k/dx;      % boundary conductance, wall is dx/2 from cell 1 centre

p = (1:nx-1)'; q = (2:nx)';

T = Tinf * ones(nx, 1);

nCheck = 6;
tCheck = linspace(tEnd/nCheck, tEnd, nCheck);
iCheck = round(tCheck/dt);
sNum  = nan(nCheck, 1);
tNum  = nan(nCheck, 1);
ci = 1;

for step = 1:nSteps
    t = step*dt;

    Tguess = T;
    for picard = 1:3
        % Energy-conserving ("secant") apparent heat capacity -- see
        % apparentCapacity.m. The plain tangent-derivative form
        % under-counts latent heat whenever a time step's temperature
        % change is comparable to the mushy-zone width, which biases
        % the melt front fast; the secant form does not.
        Capp = apparentCapacity(Tguess, T, Ts, Tl, rho*cp, rho*L);

        rows = [p; q; p; q; (1:nx)'; 1];
        cols = [q; p; p; q; (1:nx)'; 1];
        vals = [-Gc_int*ones(nx-1,1); -Gc_int*ones(nx-1,1); ...
                  Gc_int*ones(nx-1,1);  Gc_int*ones(nx-1,1); ...
                  Capp*(dx/dt); Gc_bc];
        A = sparse(rows, cols, vals, nx, nx);

        b = Capp .* (dx/dt) .* T;
        b(1) = b(1) + Gc_bc*T0;

        Tnew = A \ b;
        if norm(Tnew-Tguess, Inf) < 1e-4
            Tguess = Tnew;
            break;
        end
        Tguess = Tnew;
    end
    T = Tguess;

    if ci <= nCheck && step == iCheck(ci)
        sNum(ci) = frontPosition(xc, T, Tm);
        tNum(ci) = t;
        ci = ci + 1;
    end
end

%% ---- Compare against analytical solution -----------------------------
sAnalytic = 2*lambda*sqrt(alpha*tNum);
relErr = abs(sNum - sAnalytic) ./ sAnalytic * 100;

fprintf('\n%8s %12s %12s %10s\n', 't [s]', 's_num [m]', 's_exact [m]', 'err [%]');
for i = 1:nCheck
    fprintf('%8.1f %12.5f %12.5f %10.2f\n', tNum(i), sNum(i), sAnalytic(i), relErr(i));
end
fprintf('\nMax relative error over checkpoints: %.2f %%\n', max(relErr));
if max(relErr) < 5
    fprintf('PASS: within the 3-5%% target validation band.\n');
else
    fprintf(['WARNING: error exceeds 5%% -- refine dx/dt or shrink ', ...
             'the mushy interval dTmush.\n']);
end

%% ---- Plot --------------------------------------------------------------
outDir = fullfile(pwd, 'results');
if ~exist(outDir, 'dir'); mkdir(outDir); end

tPlot = linspace(0, tEnd, 200);
sPlot = 2*lambda*sqrt(alpha*tPlot);

f = figure('Visible', 'off');
plot(tPlot, sPlot, '-', 'LineWidth', 1.8); hold on;
plot(tNum, sNum, 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
xlabel('time [s]'); ylabel('melt front position s(t) [m]');
legend('analytical (Stefan)', 'numerical (FV, apparent C_p)', 'Location', 'northwest');
title('1-D Stefan problem validation'); grid on;
print(f, fullfile(outDir, 'fig0_stefan_validation.png'), '-dpng', '-r150');
close(f);
fprintf('\nPlot written to %s\n', fullfile(outDir, 'fig0_stefan_validation.png'));

% Note: the melt-front interpolation used inside the time loop above
% lives in the standalone helper frontPosition.m (kept as a separate
% file, not a local function, for MATLAB/Octave compatibility).
