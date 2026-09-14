function makePlotsFastMelt(P, R, Rbare, rowsNS, rowsS, chosen, band, outDir)
%MAKEPLOTSFASTMELT  Result plots for mainFastMelt.m -- mirrors
%PCM_python/fast_melt_design.py's make_plots().

% ---- fig1: gap-width sweeps (no-sleeve time, with-sleeve power) --------
f = figure('Visible', 'off', 'Position', [100 100 1000 420]);

subplot(1,2,1);
gapsNS = [rowsNS.gapMm];
tNS = [rowsNS.t999];
plot(gapsNS, tNS, 'o-', 'Color', [0.27 0.51 0.71], 'LineWidth', 1.4); hold on;
plot(xlim(), [180 180], 'r--', 'LineWidth', 1);
xlabel('PCM annulus gap width [mm]'); ylabel('full-melt time [s]');
title('Pure pipe conduction (no sleeve)');
legend('t999', '3 min target', 'Location', 'northwest');
grid on;

subplot(1,2,2);
gapsS = [rowsS.gapMm];
pS = [rowsS.Pmin];
plot(gapsS, pS, 'o-', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.4); hold on;
plot([chosen.gapMm chosen.gapMm], ylim(), 'g:', 'LineWidth', 1.5);
xlabel('PCM annulus gap width [mm]');
ylabel('min. sleeve power to melt in <170 s [W/m]');
title('With pulsed copper sleeve');
legend('P_{min}', 'recommended', 'Location', 'northwest');
grid on;

print(f, fullfile(outDir, 'fig1_gap_sweeps.png'), '-dpng', '-r150');
close(f);

% ---- fig2: liquid fraction, recommended design vs. no-sleeve at same gap --
f = figure('Visible', 'off');
plot(R.t, R.favg*100, 'LineWidth', 1.8); hold on;
plot(Rbare.t, Rbare.favg*100, '--', 'LineWidth', 1.8);
plot([180 180], ylim(), 'r:', 'LineWidth', 1);
xlabel('time [s]'); ylabel('PCM average liquid fraction [%]');
title(sprintf('Melting progress at the recommended %d mm gap', chosen.gapMm));
legend('with pulsed sleeve', 'no sleeve (same gap)', '3 min target', 'Location', 'southeast');
grid on;
print(f, fullfile(outDir, 'fig2_liquid_fraction.png'), '-dpng', '-r150');
close(f);

% ---- fig3: temperature profile snapshots ---------------------------------
G = createGeometryRadial(P, chosen.frac);
snaps = R.snapshots;
nSnap = numel(snaps);
pick = round(linspace(1, nSnap, min(7, nSnap)));
f = figure('Visible', 'off'); hold on;
cmap = jet(numel(pick));
legendEntries = {};
for i = 1:numel(pick)
    s = snaps(pick(i));
    plot(G.rc*1000, s.T, 'Color', cmap(i,:), 'LineWidth', 1.4);
    legendEntries{end+1} = sprintf('t=%.0f s', s.t); %#ok<AGROW>
end
plot(xlim(), [P.Tm P.Tm], 'r--', 'LineWidth', 1);
legendEntries{end+1} = 'PCM melt point';
yl = ylim();
patch([band(1) band(2) band(2) band(1)]*1000, [yl(1) yl(1) yl(2) yl(2)], ...
      [1 0.65 0], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
legendEntries{end+1} = 'copper sleeve';
xlabel('radius [mm]'); ylabel('temperature [\circC]');
title('Radial temperature profiles (recommended fast-melt design)');
legend(legendEntries, 'Location', 'northeast', 'FontSize', 7);
grid on;
print(f, fullfile(outDir, 'fig3_temperature_profiles.png'), '-dpng', '-r150');
close(f);

% ---- fig4: sleeve power/temperature trace --------------------------------
f = figure('Visible', 'off');
subplot(2,1,1);
plot(R.t, R.Tsleeve, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.4); hold on;
plot(xlim(), [P.Tm P.Tm], 'r--', 'LineWidth', 0.8);
ylabel('sleeve avg. temp [\circC]');
title('Recommended design: sleeve pulse-heating trace');
legend('sleeve temp', 'PCM melt point', 'Location', 'southeast');
grid on;
subplot(2,1,2);
stairs(R.t, R.power, 'Color', [0.27 0.51 0.71], 'LineWidth', 1.2);
xlabel('time [s]'); ylabel('sleeve power [W/m]');
grid on;
print(f, fullfile(outDir, 'fig4_sleeve_pulsing.png'), '-dpng', '-r150');
close(f);

end
