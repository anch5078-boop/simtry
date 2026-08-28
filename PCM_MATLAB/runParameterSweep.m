function T = runParameterSweep(r1Rratios, r2Rratios, phiList, baseP, verbose)
%RUNPARAMETERSWEEP  Batch-run the Version-1 solver over heater-radius
%ratios r1/R, r2/R and Cu-particle volume fraction phi (Section 13 of
%the design notes).
%
%   T = RUNPARAMETERSWEEP(r1Rratios, r2Rratios, phiList, baseP, verbose)
%
%       r1Rratios : vector of r1/R values, e.g. [0.20 0.30 0.40]
%       r2Rratios : vector of r2/R values, e.g. [0.55 0.70 0.85]
%       phiList   : vector of Cu-particle volume fractions, e.g.
%                   [0 0.005 0.01 0.02]
%       baseP     : (optional) base parameter struct from parameters();
%                   every combination overrides r1i/r2i/phi on top of
%                   this. Defaults to parameters().
%       verbose   : (optional, default false) print progress per run
%
%   Returns a MATLAB table T with one row per (r1/R, r2/R, phi)
%   combination and columns: r1R, r2R, phi, t50, t90, t99, Theater_max.
%   Cases with r2i <= r1o (heater shells would overlap/touch) are
%   skipped and reported with NaN results.
%
%   Example:
%       T = runParameterSweep([0.20 0.30 0.40], [0.55 0.70 0.85], ...
%                              [0 0.005 0.01 0.02]);
%       writetable(T, 'results/sweep.csv');

if nargin < 4 || isempty(baseP)
    baseP = parameters();
end
if nargin < 5
    verbose = false;
end

nR1 = numel(r1Rratios);
nR2 = numel(r2Rratios);
nPhi = numel(phiList);
nCases = nR1 * nR2 * nPhi;

r1R = zeros(nCases,1); r2R = zeros(nCases,1); phiCol = zeros(nCases,1);
t50 = nan(nCases,1); t90 = nan(nCases,1); t99 = nan(nCases,1);
TheaterMax = nan(nCases,1);
caseId = strings(nCases,1);

row = 0;
for i = 1:nR1
    for j = 1:nR2
        for m = 1:nPhi
            row = row + 1;
            r1R(row) = r1Rratios(i);
            r2R(row) = r2Rratios(j);
            phiCol(row) = phiList(m);
            caseId(row) = sprintf('r1R%.2f_r2R%.2f_phi%.3f', r1Rratios(i), r2Rratios(j), phiList(m));

            P = baseP;
            P.r1i = r1Rratios(i) * P.R;
            P.r1o = P.r1i + P.th1;
            P.r2i = r2Rratios(j) * P.R;
            P.r2o = P.r2i + P.th2;
            P.phi = phiList(m);

            if P.r2i <= P.r1o
                if verbose
                    fprintf('[skip] %s : heater shells overlap (r2i <= r1o)\n', caseId(row));
                end
                continue;
            end

            try
                G = createGeometry(P);
                C = compositeProperties(P);
                R = runSimulation(P, G, C, verbose);
                t50(row) = R.t50;
                t90(row) = R.t90;
                t99(row) = R.t99;
                TheaterMax(row) = R.TheaterMax(end);
            catch ME
                warning('runParameterSweep:caseFailed', '%s failed: %s', caseId(row), ME.message);
            end

            fprintf('[%3d/%3d] %-28s  t50=%8s  t90=%8s  t99=%8s\n', row, nCases, ...
                caseId(row), formatTimeOrNA(t50(row)), formatTimeOrNA(t90(row)), formatTimeOrNA(t99(row)));
        end
    end
end

T = table(caseId, r1R, r2R, phiCol, t50, t90, t99, TheaterMax, ...
    'VariableNames', {'case','r1_over_R','r2_over_R','phi','t50','t90','t99','Theater_max'});

end
