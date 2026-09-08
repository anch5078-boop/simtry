function m = milestonesFromSeries(t, frac)
%MILESTONESFROMSERIES  Linear-interpolated t50/t90/t95/t99 -- the times
%at which `frac` (a melt-fraction-like series, assumed non-decreasing
%enough to cross each threshold once) first reaches 50%/90%/95%/99%.
%
%   m = MILESTONESFROMSERIES(t, frac)
%
%   Returns a struct with fields t50, t90, t95, t99, each either a
%   scalar time [s] or NaN if that threshold was never reached within
%   the series (mirrors ccm_slip_model.py's None-if-not-reached
%   convention -- use isnan() to check, since a MATLAB struct field
%   cannot hold an absent value the way a Python dict value can be
%   None).

labels  = {'t50', 't90', 't95', 't99'};
threshs = [0.50,  0.90,  0.95,  0.99];

m = struct();
for k = 1:numel(labels)
    thresh = threshs(k);
    idx = find(frac >= thresh, 1, 'first');
    if isempty(idx) || idx <= 1
        m.(labels{k}) = NaN;
        continue;
    end
    f0 = frac(idx-1); f1 = frac(idx);
    t0 = t(idx-1);    t1 = t(idx);
    if f1 == f0
        fracInterp = 0;
    else
        fracInterp = (thresh - f0) / (f1 - f0);
    end
    m.(labels{k}) = t0 + fracInterp*(t1 - t0);
end

end
