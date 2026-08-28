function markMilestone(tval, yval, lbl)
%MARKMILESTONE  Plot a filled marker + text label at (tval,yval) on the
%current axes, e.g. to flag t50/t90/t99 on a liquid-fraction-vs-time
%plot. Does nothing if tval is NaN (milestone not reached).
if ~isnan(tval)
    plot(tval, yval, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
    text(tval, yval, ['  ' lbl], 'VerticalAlignment', 'bottom');
end
end
