function plotResults(R, G, P, outDir)
%PLOTRESULTS  Standard Version-1 result plots:
%   (1) average PCM liquid fraction vs time, with t50/t90/t99 marked
%   (2) max heater temperature and max domain temperature vs time
%   (3) mid-height temperature slice at the final stored snapshot
%   (4) mid-height liquid-fraction slice at the final snapshot
%
%   PLOTRESULTS(R, G, P, outDir) also saves each figure as a PNG into
%   outDir (default: current directory) so they can be inspected
%   without a display (headless / Octave-cli friendly).

if nargin < 4
    outDir = pwd;
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

lastSnap = R.snapshots(end);
kmid = max(1, round(G.nz/2));

%% ---- Figure 1: liquid fraction history ----
f1 = figure('Visible', 'off');
plot(R.t, R.favg, 'LineWidth', 1.8); hold on;
drawGuideline(0.50); drawGuideline(0.90); drawGuideline(0.99);
markMilestone(R.t50, 0.50, 't_{50}');
markMilestone(R.t90, 0.90, 't_{90}');
markMilestone(R.t99, 0.99, 't_{99}');
xlabel('time [s]'); ylabel('average PCM liquid fraction  f_l');
title('PCM melting progress');
ylim([0 1]); grid on;
print(f1, fullfile(outDir, 'fig1_liquidFraction_vs_time.png'), '-dpng', '-r150');
close(f1);

%% ---- Figure 2: temperature history ----
f2 = figure('Visible', 'off');
plot(R.t, R.TheaterMax, 'LineWidth', 1.8); hold on;
plot(R.t, R.Tmax, '--', 'LineWidth', 1.4);
xlabel('time [s]'); ylabel('temperature [^{\circ}C]');
legend('T_{heater,max}', 'T_{max} (domain)', 'Location', 'southeast');
title('Peak temperatures'); grid on;
print(f2, fullfile(outDir, 'fig2_temperature_vs_time.png'), '-dpng', '-r150');
close(f2);

%% ---- Figure 3: temperature slice, mid-height, final snapshot ----
Tfull = expandToFullGrid(lastSnap.T, G, NaN);
Tslice = squeeze(Tfull(:, :, kmid))';   % transpose -> rows=y, cols=x for plotting

f3 = figure('Visible', 'off');
im3 = imagesc(G.xc, G.yc, Tslice); axis equal tight; set(gca, 'YDir', 'normal');
set(im3, 'AlphaData', double(~isnan(Tslice)));   % render inactive (outside-cylinder) cells transparent
colorbar; xlabel('x [m]'); ylabel('y [m]');
title(sprintf('Temperature slice z \\approx %.4f m, t = %.1f s [^{\\circ}C]', ...
    G.zc(kmid), lastSnap.t));
print(f3, fullfile(outDir, 'fig3_temperature_slice.png'), '-dpng', '-r150');
close(f3);

%% ---- Figure 4: liquid fraction slice, mid-height, final snapshot ----
flFull = expandToFullGrid(lastSnap.fl, G, NaN);
flSlice = squeeze(flFull(:, :, kmid))';

f4 = figure('Visible', 'off');
im4 = imagesc(G.xc, G.yc, flSlice, [0 1]); axis equal tight; set(gca, 'YDir', 'normal');
set(im4, 'AlphaData', double(~isnan(flSlice)));   % render inactive (outside-cylinder) cells transparent
colorbar; xlabel('x [m]'); ylabel('y [m]');
title(sprintf('Liquid fraction slice z \\approx %.4f m, t = %.1f s', ...
    G.zc(kmid), lastSnap.t));
print(f4, fullfile(outDir, 'fig4_liquidFraction_slice.png'), '-dpng', '-r150');
close(f4);

fprintf('Plots written to %s\n', outDir);

end
