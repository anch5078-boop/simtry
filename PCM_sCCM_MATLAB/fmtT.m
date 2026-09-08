function s = fmtT(x)
%FMTT  Format a milestone time (seconds, or NaN if never reached) for
%the runAll.m console table.
if isnan(x)
    s = 'not reached';
else
    s = sprintf('%6.1f s (%5.2f min)', x, x/60);
end
end
