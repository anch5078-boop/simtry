%RUNALL  Driver script: runs Model A (ccmSlipModel, 0-D lumped sCCM
%force-balance) and Model B (radialEnthalpyModel, 1-D radial transient
%enthalpy PDE) across the three named slip-length scenarios grounded on
%Li et al. 2026 (Nature), prints a milestone summary table, and saves
%comparison plots to results/.
%
%Usage (MATLAB or GNU Octave, tested on Octave 8.4):
%   cd PCM_sCCM_MATLAB
%   runAll

close all;
if exist('graphics_toolkit', 'file')
    % Octave only (MATLAB has no such function). Prefer qt if available
    % -- gnuplot is only a fallback for headless environments without
    % it (see README.md's "Octave/gnuplot rendering notes").
    try
        if ~any(strcmp(available_graphics_toolkits(), 'qt'))
            graphics_toolkit('gnuplot');
        end
    catch
        % Ignore -- worst case we get whatever Octave's default is.
    end
end

thisDir = fileparts(mfilename('fullpath'));
resultsDir = fullfile(thisDir, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

P = parameters();

scenarioNames = {'no slip (b=0, baseline CCM)', ...
                  'slip b=45 um (paper, low end)', ...
                  'slip b=90 um (paper, high end)'};
scenarioB = [0.0, 45.0e-6, 90.0e-6];
colors = [0.298 0.431 0.961;    % blue
          0.961 0.624 0.000;    % orange
          0.878 0.192 0.192];   % red

fprintf('%s\n', repmat('=', 1, 78));
fprintf('sCCM vertical-tube pulsed-heater model -- Model A (0-D lumped CCM)\n');
fprintf('%s\n', repmat('=', 1, 78));

resultsA = cell(1, numel(scenarioNames));
for k = 1:numel(scenarioNames)
    resultsA{k} = ccmSlipModel(P, scenarioB(k));
end

printMilestoneTable(scenarioNames, resultsA);

baseT90 = resultsA{1}.milestones.t90;
fprintf('\n');
for k = 1:numel(scenarioNames)
    t90 = resultsA{k}.milestones.t90;
    if ~isnan(baseT90) && ~isnan(t90)
        fprintf('  %s: t90 speed-up vs no-slip = %5.2fx\n', scenarioNames{k}, baseT90/t90);
    end
end

fprintf('\n%s\n', repmat('=', 1, 78));
fprintf('sCCM vertical-tube pulsed-heater model -- Model B (1-D radial enthalpy PDE)\n');
fprintf('%s\n', repmat('=', 1, 78));

resultsB = cell(1, numel(scenarioNames));
for k = 1:numel(scenarioNames)
    resultsB{k} = radialEnthalpyModel(P, scenarioB(k));
end

printMilestoneTable(scenarioNames, resultsB);

%% ---- Plot 1: Model A melt fraction vs time --------------------------
fig = figure('visible', 'off');
hold on;
for k = 1:numel(scenarioNames)
    R = resultsA{k};
    plot(R.t/60, 100*R.meltFraction, 'Color', colors(k,:), 'LineWidth', 1.8, ...
         'DisplayName', scenarioNames{k});
end
xlabel('time [min]');
ylabel('PCM melted [%] (Model A, sinking-column mass basis)');
title('Model A: sCCM lumped force-balance -- melt progress');
legend(scenarioNames, 'Location', 'southeast', 'FontSize', 8);
grid on;
print(fig, fullfile(resultsDir, 'A_melt_fraction_vs_time.png'), '-dpng', '-r150');
close(fig);

%% ---- Plot 2: Model A film thickness & heat flux (first 10 min) -----
fig = figure('visible', 'off');
tcut = 600;
subplot(2,1,1); hold on;
for k = 1:numel(scenarioNames)
    R = resultsA{k};
    mask = R.t <= tcut;
    plot(R.t(mask), R.delta(mask)*1e6, 'Color', colors(k,:), 'LineWidth', 1.4, ...
         'DisplayName', scenarioNames{k});
end
ylabel('film thickness delta [um]');
title('Model A: CCM film thickness and heater-surface heat flux (first 10 min)');
legend(scenarioNames, 'Location', 'northeast', 'FontSize', 8);
grid on;

subplot(2,1,2); hold on;
base = resultsA{1};
mask = base.t <= tcut;
off = ~base.pulseOn(mask);
tOff = base.t(mask);
yMaxFlux = 0;
for k = 1:numel(scenarioNames)
    R = resultsA{k};
    mask2 = R.t <= tcut;
    yMaxFlux = max(yMaxFlux, max(R.qH(mask2)/1000));
end
shadeIntervals(tOff, off, [0 yMaxFlux*1.05]);
for k = 1:numel(scenarioNames)
    R = resultsA{k};
    mask2 = R.t <= tcut;
    % stairs(), not plot(): qH switches instantaneously between 0 (OFF)
    % and its ON value, so a step plot renders the true square wave
    % instead of a straight-line-interpolated ramp between samples.
    stairs(R.t(mask2), R.qH(mask2)/1000, 'Color', colors(k,:), 'LineWidth', 1.4);
end
ylabel('heat flux q_h [kW/m^2]');
xlabel('time [s]');
grid on;
print(fig, fullfile(resultsDir, 'A_film_and_flux.png'), '-dpng', '-r150');
close(fig);

%% ---- Plot 3: Model B liquid fraction vs time ------------------------
fig = figure('visible', 'off');
hold on;
for k = 1:numel(scenarioNames)
    R = resultsB{k};
    plot(R.t/60, 100*R.liquidFracVol, 'Color', colors(k,:), 'LineWidth', 1.8, ...
         'DisplayName', scenarioNames{k});
end
xlabel('time [min]');
ylabel('volume-averaged liquid fraction [%] (Model B, fixed annulus)');
title('Model B: 1-D radial enthalpy PDE -- melt progress');
legend(scenarioNames, 'Location', 'southeast', 'FontSize', 8);
grid on;
print(fig, fullfile(resultsDir, 'B_liquid_fraction_vs_time.png'), '-dpng', '-r150');
close(fig);

%% ---- Plot 4: Model B melt front vs time, pulsing visible ------------
fig = figure('visible', 'off');
hold on;
tcutB = 300;
baseB = resultsB{1};
maskB = baseB.t <= tcutB;
offB = ~baseB.pulseOn(maskB);
yMaxGap = 0;
for k = 1:numel(scenarioNames)
    R = resultsB{k};
    mask2 = R.t <= tcutB;
    yMaxGap = max(yMaxGap, max((R.rMelt(mask2) - P.rH)*1000));
end
shadeIntervals(baseB.t(maskB), offB, [0 yMaxGap*1.1]);
for k = 1:numel(scenarioNames)
    R = resultsB{k};
    mask2 = R.t <= tcutB;
    plot(R.t(mask2), (R.rMelt(mask2) - P.rH)*1000, 'Color', colors(k,:), ...
         'LineWidth', 1.6, 'DisplayName', scenarioNames{k});
end
xlabel('time [s]');
ylabel('melt-layer thickness  r_{melt} - r_H  [mm]');
title('Model B: melt-front growth, pulsing stalls visible (first 5 min)');
legend(scenarioNames, 'Location', 'northwest', 'FontSize', 8);
grid on;
print(fig, fullfile(resultsDir, 'B_melt_front_pulsing.png'), '-dpng', '-r150');
close(fig);

%% ---- Plot 5: Model B temperature-field snapshots (b=45 um) ----------
fig = figure('visible', 'off', 'Position', [100 100 950 500]);
hold on;
R = resultsB{2};   % the paper's low-end slip scenario, b=45 um
nSnap = numel(R.THistory);
idxs = round(linspace(1, nSnap, min(8, nSnap)));
% Simple dark-blue -> dark-red gradient (avoids jet()'s washed-out pale
% band in the middle of the range, which is hard to read on white).
nc = numel(idxs);
cmap = [linspace(0.10, 0.75, nc)', linspace(0.10, 0.10, nc)', linspace(0.75, 0.10, nc)'];
snapLabels = cell(1, numel(idxs));
snapHandles = zeros(1, numel(idxs));
for j = 1:numel(idxs)
    idx = idxs(j);
    snapHandles(j) = plot((R.rC - P.rH)*1000, R.THistory{idx}, 'Color', cmap(j,:));
    snapLabels{j} = sprintf('t=%.0f s', R.tSaved(idx));
end
hTs = plot(xlim(), [P.Ts P.Ts], 'k:', 'LineWidth', 0.8);
plot(xlim(), [P.Tl P.Tl], 'k:', 'LineWidth', 0.8);
xlabel('radial distance from heater surface  r - r_H  [mm]');
ylabel('temperature [deg C]');
title('Model B: radial temperature profiles (slip b=45 um (paper, low end))');
legend([snapHandles, hTs], [snapLabels, {'Ts, Tl (mushy band)'}], ...
       'Location', 'eastoutside', 'FontSize', 7);
grid on;
print(fig, fullfile(resultsDir, 'B_temperature_profiles.png'), '-dpng', '-r150');
close(fig);

fprintf('\nPlots written to %s\n', resultsDir);
