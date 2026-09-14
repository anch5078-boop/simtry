function s = formatTimeOrNARadial(x)
%FORMATTIMEORNARADIAL  "Hh Mm Ss (N s)" for a finite scalar, else "n/a".
if isempty(x) || isnan(x)
    s = 'n/a';
    return;
end
h = floor(x/3600);
m = floor(mod(x,3600)/60);
sec = mod(x,60);
s = sprintf('%dh %02dm %02ds (%.0f s)', h, m, sec, x);
end
