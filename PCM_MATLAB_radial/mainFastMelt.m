%MAINFASTMELT  Re-dimension the pipe/PCM/sleeve/wall annulus so the
%whole PCM charge melts in under 3 minutes, keeping the pipe radius and
%PCM material fixed. MATLAB/Octave counterpart of
%PCM_python/fast_melt_design.py -- see PCM_MATLAB_radial/README.md for
%the full write-up.
%
%   Run directly (MATLAB or Octave):  >> mainFastMelt

clear; clc; close all;
addpath(pwd);

outDir = fullfile(pwd, 'results_fast_melt');
if ~exist(outDir, 'dir'); mkdir(outDir); end

TARGET_S = 180.0;         % the <3-minute requirement
BISECT_TARGET_S = 170.0;  % bisect inside the target for margin
T_SAFETY_CAP = 180.0;     % deg C, sleeve heater cuts off above this
REALISTIC_POWER_CAP = 2000.0;  % W/m -- "realistic" cutoff for the recommended design

Pbase = parametersRadial();
rPipe = Pbase.rPipe;

fprintf('=== Fast-melt design: re-dimension for <3-minute melting ===\n\n');

%% ---- Step 1: pure pipe conduction (no sleeve) --------------------------
fprintf('--- Step 1: pure pipe conduction (no sleeve) -- gap-width sweep ---\n');
gapsNoSleeve = [2 3 4 5 6 8 10 12];
rowsNS = sweepNoSleeveGaps(rPipe, gapsNoSleeve, 220, 0.25, 1200, Pbase, true);

%% ---- Step 2: with pulsed copper sleeve ----------------------------------
fprintf('\n--- Step 2: with pulsed copper sleeve -- gap-width sweep (target <=%.0f s, cap %.0f C) ---\n', ...
    BISECT_TARGET_S, T_SAFETY_CAP);
gapsSleeve = [6 7 8 9 10 11 12 13 14];
rowsS = sweepWithSleeveGaps(rPipe, gapsSleeve, 260, 0.2, 600, 0.001, Pbase, ...
    BISECT_TARGET_S, T_SAFETY_CAP, true);

%% ---- Pick the recommended design ----------------------------------------
feasible = rowsS(~isnan([rowsS.Pmin]));
candidates = feasible([feasible.Pmin] <= REALISTIC_POWER_CAP);
[~, iMax] = max([candidates.gapMm]);
chosen = candidates(iMax);
fprintf('\nRecommended design: gap=%d mm, frac=%.2f, Pmin=%.1f W/m -> t999=%.1f s\n', ...
    chosen.gapMm, chosen.frac, chosen.Pmin, chosen.t999);

gapMm = chosen.gapMm;
frac = chosen.frac;
Pdesign = round((chosen.Pmin * 1.05)/10)*10;  % +5% margin, rounded to 10 W/m

P = Pbase;
P.rPipe = rPipe; P.rWall = rPipe + gapMm/1000;
P.nCells = 400; P.dt = 0.1; P.tMax = 400; P.tSleeve = 0.001;
matp = materialStruct(P);

G = createGeometryRadial(P, frac);
R = runFullMeltRadial(P, G, matp, Pdesign, T_SAFETY_CAP, 20, true);
fprintf('Final verified run (finer grid, +5%% power margin = %.0f W/m): t50=%.1fs t90=%.1fs t99=%.1fs t999=%.1fs\n', ...
    Pdesign, R.t50, R.t90, R.t99, R.t999);

%% ---- No-sleeve check at the same gap, for contrast ------------------------
Gbare = createGeometryNoSleeve(P);
Rbare = runFullMeltRadial(P, Gbare, matp, 0, P.Tm, 20, false);

%% ---- Save CSVs --------------------------------------------------------------
fid = fopen(fullfile(outDir, 'no_sleeve_sweep.csv'), 'w');
fprintf(fid, 'gap_mm,t999\n');
for i = 1:numel(rowsNS)
    fprintf(fid, '%d,%s\n', rowsNS(i).gapMm, numOrEmpty(rowsNS(i).t999));
end
fclose(fid);

fid = fopen(fullfile(outDir, 'sleeve_sweep.csv'), 'w');
fprintf(fid, 'gap_mm,frac,P_min_W_per_m,t999\n');
for i = 1:numel(rowsS)
    fprintf(fid, '%d,%.4f,%s,%s\n', rowsS(i).gapMm, rowsS(i).frac, ...
        numOrEmpty(rowsS(i).Pmin), numOrEmpty(rowsS(i).t999));
end
fclose(fid);

fid = fopen(fullfile(outDir, 'summary.txt'), 'w');
fprintf(fid, 'target_s=%.1f\n\n[recommended]\nr_pipe_mm=%.4f\ngap_mm=%d\nr_wall_mm=%.4f\n', ...
    TARGET_S, P.rPipe*1000, gapMm, P.rWall*1000);
fprintf(fid, 'sleeve_thickness_mm=%.4f\nsleeve_center_frac=%.4f\nsleeve_band_mm=%.4f-%.4f\n', ...
    P.tSleeve*1000, frac, G.rSleeveIn*1000, G.rSleeveOut*1000);
fprintf(fid, 'P_design_W_per_m=%.1f\nT_safety_cap_C=%.1f\n', Pdesign, T_SAFETY_CAP);
fprintf(fid, 't50=%.2f\nt90=%.2f\nt99=%.2f\nt999=%.2f\n', R.t50, R.t90, R.t99, R.t999);
fclose(fid);

%% ---- Plots --------------------------------------------------------------
makePlotsFastMelt(P, R, Rbare, rowsNS, rowsS, chosen, [G.rSleeveIn, G.rSleeveOut], outDir);

fprintf('\nWrote plots, CSVs and summary.txt to %s/\n', outDir);
fprintf('\n=== FINAL RECOMMENDATION ===\n');
fprintf('Pipe radius (fixed, given): %.2f mm\n', P.rPipe*1000);
fprintf('Container inner radius: %.2f mm (PCM annulus thickness: %d mm, was 32.3 mm)\n', ...
    P.rWall*1000, gapMm);
fprintf('Copper sleeve: %.1f mm thick, centered at %.0f%% of the way from pipe to wall (%.2f-%.2f mm)\n', ...
    P.tSleeve*1000, frac*100, G.rSleeveIn*1000, G.rSleeveOut*1000);
fprintf('Sleeve pulse power: %.0f W/m of pipe length (e.g. %.0f W for a 0.3 m section), capped off above %.0f C\n', ...
    Pdesign, Pdesign*0.3, T_SAFETY_CAP);
fprintf('Result: fully melted (99.9%%) in %.0f s = %.2f min (target: < 180 s)\n', R.t999, R.t999/60);
fprintf('Without the sleeve at this same gap (400 s cutoff): t50=%s t90=%s t99=%s t999=%s\n', ...
    numOrEmpty(Rbare.t50), numOrEmpty(Rbare.t90), numOrEmpty(Rbare.t99), numOrEmpty(Rbare.t999));
