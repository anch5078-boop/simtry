function makePlotsAsDrawn(P, Rbase, Rref, Ropt, rowsSorted, best, bandBase, bandOpt, outDir)
%MAKEPLOTSASDRAWN  Result plots for mainAsDrawn.m -- mirrors
%PCM_python/main.py's make_plots(): melting-progress comparison,
%position-sweep curve, temperature-profile snapshots, sleeve pulsing
%trace.

% ---- fig1: liquid fraction history -----------------------------------
f = figure('Visible', 'off');
plot(Rbase.t/60, Rbase.favg*100, '-', 'LineWidth', 1.8); hold on;
plot(Ropt.t/60,  Ropt.favg*100,  '-', 'LineWidth', 1.8);
plot(Rref.t/60,  Rref.favg*100,  '--', 'LineWidth', 1.8);
yline_manual(99);
xlabel('time [min]'); ylabel('PCM average liquid fraction [%]');
title('Melting progress: effect of the pulsed copper sleeve');
legend('Sleeve centered in gap (as drawn)', ...
       sprintf('Optimal position (%.0f%% pipe->wall)', best.frac*100), ...
       'No sleeve (pipe conduction only)', 'Location', 'southeast');
grid on;
print(f, fullfile(outDir, 'fig1_liquid_fraction.png'), '-dpng', '-r150');
close(f);

% ---- fig2: sweep curve -------------------------------------------------
fracs = [rowsSorted.frac];
t99s  = [rowsSorted.t99]/60;
t999s = [rowsSorted.t999]/60;
f = figure('Visible', 'off');
plot(fracs, t99s, 'o-', 'MarkerSize', 4); hold on;
plot(fracs, t999s, 's-', 'MarkerSize', 4);
xline_manual(best.frac, 'g');
xline_manual(0.5, [0.5 0.5 0.5]);
xlabel('sleeve position: fraction of pipe->wall gap');
ylabel('melting time [min]');
title('Melting time vs. copper-sleeve radial position');
legend('t99 (99% melted)', 't99.9 (fully melted)', 'optimum', 'centered (as drawn)', ...
       'Location', 'northeast');
grid on;
print(f, fullfile(outDir, 'fig2_position_sweep.png'), '-dpng', '-r150');
close(f);

% ---- fig3: temperature profile snapshots -------------------------------
G = createGeometryRadial(P, best.frac);
snaps = Ropt.snapshots;
nSnap = numel(snaps);
pick = round(linspace(1, nSnap, min(7, nSnap)));
f = figure('Visible', 'off'); hold on;
cmap = jet(numel(pick));
legendEntries = {};
for i = 1:numel(pick)
    s = snaps(pick(i));
    plot(G.rc*1000, s.T, 'Color', cmap(i,:), 'LineWidth', 1.4);
    legendEntries{end+1} = sprintf('t=%.0f min', s.t/60); %#ok<AGROW>
end
plot(xlim(), [P.Tm P.Tm], 'r--', 'LineWidth', 1);
legendEntries{end+1} = 'PCM melt point';
patch([bandOpt(1) bandOpt(2) bandOpt(2) bandOpt(1)]*1000, ...
      [min(ylim()) min(ylim()) max(ylim()) max(ylim())], ...
      [1 0.65 0], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
legendEntries{end+1} = 'copper sleeve';
xlabel('radius [mm]'); ylabel('temperature [\circC]');
title('Radial temperature profiles (optimal sleeve position)');
legend(legendEntries, 'Location', 'northwest', 'FontSize', 7);
grid on;
print(f, fullfile(outDir, 'fig3_temperature_profiles.png'), '-dpng', '-r150');
close(f);

% ---- fig4: sleeve temperature + pulsed power (baseline, centered) -------
f = figure('Visible', 'off');
subplot(2,1,1);
plot(Rbase.t/60, Rbase.Tsleeve, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.4); hold on;
plot(xlim(), [P.Tm P.Tm], 'r--', 'LineWidth', 1);
ylabel('sleeve avg. temp [\circC]');
title('Sleeve pulse-heating behaviour (centered case)');
grid on;
subplot(2,1,2);
stairs(Rbase.t/60, Rbase.power, 'Color', [0.27 0.51 0.71], 'LineWidth', 1.2);
xlabel('time [min]'); ylabel('sleeve power [W/m]');
grid on;
print(f, fullfile(outDir, 'fig4_sleeve_pulsing.png'), '-dpng', '-r150');
close(f);

end

% ------------------------------------------------------------------------
function yline_manual(y)
xl = xlim();
plot(xl, [y y], ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.7);
end

function xline_manual(x, color)
yl = ylim();
plot([x x], yl, ':', 'Color', color, 'LineWidth', 1.2);
end
