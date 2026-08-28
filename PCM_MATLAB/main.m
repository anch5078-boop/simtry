%MAIN  Version-1 driver script: 3D transient conduction + PCM melting
%(apparent heat capacity) in a cylinder heated by two concentric
%copper heater sleeves, with uniform Cu-particle enhancement of the
%PCM and optional pulsed heating.
%
%   Run this script directly (MATLAB or Octave):
%       >> main
%
%   It builds the geometry, runs the simulation, prints t50/t90/t99,
%   and writes result plots (PNG) to ./results/.

clear; clc; close all;

addpath(pwd);   % make sure the PCM_MATLAB functions are on the path

fprintf('=== PCM Version-1 solver ===\n');

%% ---- 1) Parameters -------------------------------------------------
P = parameters();

%% ---- 2) Geometry -----------------------------------------------------
fprintf('Building geometry (h = %.4g m) ...\n', P.h);
G = createGeometry(P);
fprintf('  grid: %d x %d x %d, %d active cells (%d PCM, %d heater1, %d heater2)\n', ...
    G.nx, G.ny, G.nz, G.Na, nnz(G.pcmMaskA), nnz(G.heater1MaskA), nnz(G.heater2MaskA));

%% ---- 3) Composite (Cu-particle enhanced) PCM properties -------------
C = compositeProperties(P);
fprintf('  phi = %.4f | ks_eff = %.4f, kl_eff = %.4f W/mK | (rho L)_eff = %.1f J/m^3\n', ...
    P.phi, C.ks_eff, C.kl_eff, C.rhoL_eff);

%% ---- 4) Run the transient simulation ----------------------------------
fprintf('Running simulation (dt = %.2f s, tMax = %.1f s) ...\n', P.dt, P.tMax);
tic;
R = runSimulation(P, G, C, true);
fprintf('Elapsed wall time: %.1f s\n', toc);

%% ---- 5) Report ---------------------------------------------------------
fprintf('\n--- Results ---\n');
fprintf('t50 = %s s\n', formatTimeOrNA(R.t50));
fprintf('t90 = %s s\n', formatTimeOrNA(R.t90));
fprintf('t99 = %s s\n', formatTimeOrNA(R.t99));
fprintf('Final T_heater,max = %.2f C\n', R.TheaterMax(end));

%% ---- 6) Plots ------------------------------------------------------
plotResults(R, G, P, fullfile(pwd, 'results'));

%% ---- 7) Save raw results ------------------------------------------
save(fullfile(pwd, 'results', 'run.mat'), 'R', 'G', 'P', 'C');
fprintf('Saved results/run.mat\n');
