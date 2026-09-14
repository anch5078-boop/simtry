function s = numOrEmpty(x)
%NUMOREMPTY  "%.4f"-formatted number, or an empty string for NaN
%(kept as a separate file, not a script-local function, for
%MATLAB/Octave compatibility -- mirrors the convention already used by
%PCM_MATLAB/formatTimeOrNA.m et al.).
if isnan(x)
    s = '';
else
    s = sprintf('%.4f', x);
end
end
