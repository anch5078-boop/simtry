%MAINASDRAWN  Driver script: the case as drawn (copper sleeve centered
%in the pipe->wall gap), pulse-heated at the lowest power that reaches
%the PCM melt point within the target ramp time; a no-sleeve (pipe-
%only) reference; and a sweep of the sleeve's radial position to find
%the minimum achievable melting time for this architecture. This is
%the MATLAB/Octave counterpart of PCM_python/main.py -- see
%PCM_MATLAB_radial/README.md for the full write-up.
%
%   Run directly (MATLAB or Octave):  >> mainAsDrawn

clear; clc; close all;
addpath(pwd);

outDir = fullfile(pwd, 'results');
if ~exist(outDir, 'dir'); mkdir(outDir); end

P = parametersRadial();
matp = materialStruct(P);

fprintf('=== PCM radial model (MATLAB/Octave): hot pipe + pulsed copper sleeve ===\n');
fprintf('rPipe=%.2f mm, rWall=%.2f mm (gap %.2f mm) [container size ASSUMED]\n', ...
    P.rPipe*1000, P.rWall*1000, (P.rWall-P.rPipe)*1000);
fprintf('Thotwater=%.0f C, Tm=%.0f C, Tinit=%.0f C\n', P.Thotwater, P.Tm, P.Tinit);
fprintf('PCM: rho=%.0f, cp=%.0f J/kgK, ks=%.2f, kl=%.2f W/mK, L=%.0f J/kg [ASSUMED dataset]\n\n', ...
    P.rhoPCM, P.cpPCM, P.kPCMs, P.kPCMl, P.LPCM);

%% ---- 1) Baseline: sleeve centered in the gap, as drawn ------------------
fprintf('--- Baseline: sleeve centered in the gap (as drawn) ---\n');
Gbase = createGeometryRadial(P, 0.5);
tic;
[PminBase, TreachedBase] = findMinPowerRadial(P, Gbase, matp);
Rbase = runFullMeltRadial(P, Gbase, matp, PminBase, P.Tm, 20, true);
fprintf('  sleeve band: %.2f-%.2f mm\n', Gbase.rSleeveIn*1000, Gbase.rSleeveOut*1000);
fprintf('  minimum sleeve peak power Pmin = %.1f W/m (reaches %.1f C by tRamp=%.0f s)\n', ...
    PminBase, TreachedBase, P.tRampTarget);
fprintf('  t50=%s  t90=%s  t99=%s  t999(full melt)=%s\n', ...
    formatTimeOrNARadial(Rbase.t50), formatTimeOrNARadial(Rbase.t90), ...
    formatTimeOrNARadial(Rbase.t99), formatTimeOrNARadial(Rbase.t999));
fprintf('  sleeve heater duty cycle: %.1f%%\n', Rbase.dutyCycle*100);
fprintf('  elapsed wall time: %.1f s\n\n', toc);

%% ---- 2) No-sleeve reference ----------------------------------------------
fprintf('--- Reference: no sleeve, pipe conduction only ---\n');
Gref = createGeometryNoSleeve(P);
Rref = runFullMeltRadial(P, Gref, matp, 0, P.Tm, 20, false);
fprintf('  t50=%s  t90=%s  t99=%s  t999(full melt)=%s\n\n', ...
    formatTimeOrNARadial(Rref.t50), formatTimeOrNARadial(Rref.t90), ...
    formatTimeOrNARadial(Rref.t99), formatTimeOrNARadial(Rref.t999));

%% ---- 3) Sleeve-position sweep --------------------------------------------
fprintf('--- Sleeve-position sweep (coarse) ---\n');
fracsCoarse = linspace(0.15, 0.90, 16);
rowsCoarse = runPositionSweepRadial(P, matp, fracsCoarse, 40, true);

[~, iBest] = min([rowsCoarse.t999]);
bestCoarse = rowsCoarse(iBest);
fprintf('\n  coarse best: frac=%.2f -> t999=%.0f s\n', bestCoarse.frac, bestCoarse.t999);

fprintf('\n--- Sleeve-position sweep (refine near coarse optimum) ---\n');
lo = max(0.10, bestCoarse.frac - 0.08);
hi = min(0.95, bestCoarse.frac + 0.08);
fracsFine = linspace(lo, hi, 9);
rowsFine = runPositionSweepRadial(P, matp, fracsFine, 40, true);

rowsAll = [rowsCoarse, rowsFine];
[~, iBestAll] = min([rowsAll.t999]);
best = rowsAll(iBestAll);

fprintf('\n  OPTIMAL sleeve position: frac=%.3f (rMid=%.2f mm)\n', best.frac, best.rMid*1000);
fprintf('  Pmin at optimum = %.1f W/m\n', best.Pmin);
fprintf('  t99=%s   t999(full melt)=%s\n', formatTimeOrNARadial(best.t99), formatTimeOrNARadial(best.t999));

%% ---- Full run at the optimal position (for plotting) ----------------------
Gopt = createGeometryRadial(P, best.frac);
Ropt = runFullMeltRadial(P, Gopt, matp, best.Pmin, P.Tm, 20, true);

%% ---- Save CSV of the sweep -------------------------------------------------
[~, order] = sort([rowsAll.frac]);
rowsSorted = rowsAll(order);
fid = fopen(fullfile(outDir, 'sweep.csv'), 'w');
fprintf(fid, 'frac,rIn,rOut,rMid,Pmin_W_per_m,t50,t90,t99,t999,dutyCycle,converged\n');
for i = 1:numel(rowsSorted)
    r = rowsSorted(i);
    fprintf(fid, '%.4f,%.6f,%.6f,%.6f,%.4f,%s,%s,%s,%s,%.4f,%d\n', ...
        r.frac, r.rIn, r.rOut, r.rMid, r.Pmin, ...
        numOrEmpty(r.t50), numOrEmpty(r.t90), numOrEmpty(r.t99), numOrEmpty(r.t999), ...
        r.dutyCycle, r.converged);
end
fclose(fid);

%% ---- Plots -----------------------------------------------------------------
makePlotsAsDrawn(P, Rbase, Rref, Ropt, rowsSorted, best, ...
    [Gbase.rSleeveIn, Gbase.rSleeveOut], [Gopt.rSleeveIn, Gopt.rSleeveOut], outDir);

%% ---- Summary -----------------------------------------------------------------
fid = fopen(fullfile(outDir, 'summary.txt'), 'w');
fprintf(fid, 'rPipe_mm=%.4f\nrWall_mm=%.4f\nThotwater=%.2f\nTm=%.2f\nTinit=%.2f\ntRampTarget=%.1f\n\n', ...
    P.rPipe*1000, P.rWall*1000, P.Thotwater, P.Tm, P.Tinit, P.tRampTarget);
fprintf(fid, '[baseline_centered]\nPmin_W_per_m=%.4f\nt50=%s\nt90=%s\nt99=%s\nt999=%s\ndutyCycle=%.4f\n\n', ...
    PminBase, numOrEmpty(Rbase.t50), numOrEmpty(Rbase.t90), numOrEmpty(Rbase.t99), numOrEmpty(Rbase.t999), Rbase.dutyCycle);
fprintf(fid, '[no_sleeve_reference]\nt50=%s\nt90=%s\nt99=%s\nt999=%s\n\n', ...
    numOrEmpty(Rref.t50), numOrEmpty(Rref.t90), numOrEmpty(Rref.t99), numOrEmpty(Rref.t999));
fprintf(fid, '[optimal]\nfrac=%.4f\nrMid_mm=%.4f\nPmin_W_per_m=%.4f\nt99=%s\nt999=%s\n', ...
    best.frac, best.rMid*1000, best.Pmin, numOrEmpty(best.t99), numOrEmpty(best.t999));
fclose(fid);

fprintf('\nWrote plots, sweep.csv and summary.txt to %s/\n', outDir);
fprintf('\n=== HEADLINE RESULT ===\n');
fprintf('As drawn (sleeve centered, %.0f W/m minimum pulse power): full melt in %s\n', ...
    PminBase, formatTimeOrNARadial(Rbase.t999));
fprintf('No auxiliary sleeve heater at all (pipe conduction only): full melt in %s\n', ...
    formatTimeOrNARadial(Rref.t999));
fprintf('Minimum achievable melting time (sleeve at %.0f%% pipe->wall, %.0f W/m): full melt in %s\n', ...
    best.frac*100, best.Pmin, formatTimeOrNARadial(best.t999));
