function s = formatTimeOrNA(x)
%FORMATTIMEORNA  Small display helper: "%.2f" for a finite value, or a
%readable "not reached" label for NaN. Used by main.m and
%runParameterSweep.m when reporting t50/t90/t99.
if isnan(x)
    s = 'NaN (not reached)';
else
    s = sprintf('%.2f', x);
end
end
