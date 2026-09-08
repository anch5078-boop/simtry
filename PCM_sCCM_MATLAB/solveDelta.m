function delta = solveDelta(Wnet, muL, kL, dT, Hc, w, rhoL, Lf, b)
%SOLVEDELTA  Root-find the quasi-steady CCM film thickness delta [m]
%satisfying the lubrication force balance derived in ccmSlipModel.m:
%
%   Wnet = 4*muL*w*kL*dT*Hc^3 ./ (delta .* deltaEffCubed(delta,b) * rhoL*Lf)
%
%given the current driving dT and wetted height Hc. Solved with fzero
%bracketed on [1e-9, 5e-2] m (the same bounds ccm_slip_model.py's
%brentq call uses).
%
%   delta = SOLVEDELTA(Wnet, muL, kL, dT, Hc, w, rhoL, Lf, b)
%
%   Returns the bracket endpoint (capped) if no root exists in range --
%   this occurs only as Wnet -> 0 (near complete melt-out; see
%   ccmSlipModel.m notes) or, at the thin end, if even the thinnest
%   allowed film cannot support the current weight at this dT.

lo = 1e-9;
hi = 5e-2;

predicted = @(d) 4*muL*w*Hc^3*kL*dT ./ (d .* deltaEffCubed(d, b) .* rhoL*Lf);

fLo = predicted(lo) - Wnet;
fHi = predicted(hi) - Wnet;

if fLo <= 0
    delta = lo;
    return;
end
if fHi >= 0
    delta = hi;
    return;
end

delta = fzero(@(d) predicted(d) - Wnet, [lo, hi]);

end
