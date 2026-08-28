function drawGuideline(y)
%DRAWGUIDELINE  Draw a light dotted horizontal reference line at
%height y across the current axes' x-limits (e.g. the 0.5/0.9/0.99
%liquid-fraction thresholds).
xl = xlim();
plot(xl, [y y], ':', 'Color', [0.6 0.6 0.6]);
end
