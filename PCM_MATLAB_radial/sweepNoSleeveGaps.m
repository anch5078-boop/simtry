function rows = sweepNoSleeveGaps(rPipe, gapListMm, nCells, dt, tMax, Pbase, verbose)
%SWEEPNOSLEEVEGAPS  Pure pipe-conduction (no sleeve) full-melt time vs.
%PCM annulus gap width -- shows the "no auxiliary heater needed" limit
%for the <3-minute target.

if nargin < 7 || isempty(verbose), verbose = true; end

rows = struct('gapMm', {}, 't999', {});
for i = 1:numel(gapListMm)
    gapMm = gapListMm(i);
    P = Pbase;
    P.rPipe = rPipe;
    P.rWall = rPipe + gapMm/1000;
    P.nCells = nCells; P.dt = dt; P.tMax = tMax;
    matp = materialStruct(P);
    G = createGeometryNoSleeve(P);
    R = runFullMeltRadial(P, G, matp, 0, P.Tm, 40, false);

    row.gapMm = gapMm; row.t999 = R.t999;
    rows(end+1) = row; %#ok<AGROW>
    if verbose
        fprintf('  [no sleeve] gap=%3d mm  t999=%s\n', gapMm, numOrNAstr(R.t999));
    end
end

end

function s = numOrNAstr(x)
if isnan(x), s = 'n/a'; else, s = sprintf('%.2f', x); end
end
