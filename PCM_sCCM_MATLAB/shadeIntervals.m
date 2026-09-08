function shadeIntervals(t, offMask, yRange)
%SHADEINTERVALS  Shade contiguous True runs of offMask (both length
%numel(t)) as grey vertical bands spanning yRange, used by runAll.m to
%mark pulse-OFF intervals on a time-series plot. Call this before
%plotting the data series so the shading sits behind the lines.
if ~any(offMask)
    return;
end
d = diff([false; offMask(:); false]);
starts = find(d == 1);
ends   = find(d == -1) - 1;
for i = 1:numel(starts)
    x0 = t(starts(i));
    x1 = t(min(ends(i), numel(t)));
    fill([x0 x1 x1 x0], [yRange(1) yRange(1) yRange(2) yRange(2)], ...
         [0.85 0.85 0.85], 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
end
