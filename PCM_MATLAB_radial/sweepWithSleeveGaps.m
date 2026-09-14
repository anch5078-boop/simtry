function rows = sweepWithSleeveGaps(rPipe, gapListMm, nCells, dt, tMax, tSleeve, Pbase, targetS, Tsetpoint, verbose)
%SWEEPWITHSLEEVEGAPS  With the pulsed copper sleeve: for each PCM
%annulus gap width, re-optimize the sleeve's radial position (quick
%ranking at a fixed reference power) and bisect for the minimum sleeve
%peak power that still melts everything within targetS seconds
%(sleeve bang-bang capped at Tsetpoint for safety).

if nargin < 8  || isempty(targetS),   targetS = 170.0; end
if nargin < 9  || isempty(Tsetpoint), Tsetpoint = 180.0; end
if nargin < 10 || isempty(verbose),   verbose = true; end

rows = struct('gapMm', {}, 'frac', {}, 'Pmin', {}, 't999', {});
for i = 1:numel(gapListMm)
    gapMm = gapListMm(i);
    P = Pbase;
    P.rPipe = rPipe;
    P.rWall = rPipe + gapMm/1000;
    P.nCells = nCells; P.dt = dt; P.tMax = tMax; P.tSleeve = tSleeve;
    matp = materialStruct(P);

    frac = positionSweepQuick(P, matp, linspace(0.3, 0.85, 12), 1500.0, Tsetpoint);
    G = createGeometryRadial(P, frac);
    [Pmin, t999] = minPowerForTarget(P, G, matp, targetS, Tsetpoint);

    row.gapMm = gapMm; row.frac = frac; row.Pmin = Pmin; row.t999 = t999;
    rows(end+1) = row; %#ok<AGROW>
    if verbose
        fprintf('  [with sleeve] gap=%3d mm  frac=%.2f  Pmin=%s  t999=%s\n', ...
            gapMm, frac, numOrNAstr(Pmin), numOrNAstr(t999));
    end
end

end

function s = numOrNAstr(x)
if isnan(x), s = 'n/a'; else, s = sprintf('%.2f', x); end
end
