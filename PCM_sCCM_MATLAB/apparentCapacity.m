function Capp = apparentCapacity(Tguess, Tn, Ts, Tl, rhocp, rhoL)
%APPARENTCAPACITY  Energy-conserving ("secant") apparent heat capacity
%for the phase-change term -- identical method to
%PCM_MATLAB/apparentCapacity.m (see that file for the full derivation
%of why the naive tangent-derivative form silently under-counts latent
%heat and this secant form is required for exact energy conservation):
%
%   rhocp*(Tguess-Tn) + rhoL*(fl(Tguess)-fl(Tn)) = Capp*(Tguess-Tn)
%
%   => Capp = rhocp + rhoL * (fl(Tguess)-fl(Tn)) / (Tguess-Tn)
%
%   Capp = APPARENTCAPACITY(Tguess, Tn, Ts, Tl, rhocp, rhoL)
%
%   All of Tguess, Tn may be arrays (same size); rhocp, rhoL are
%   scalars (or arrays the same size). Falls back to the tangent
%   derivative dfl/dT where Tguess and Tn coincide to machine
%   precision (secant undefined).

dT = Tguess - Tn;
tinyStep = abs(dT) < 1e-8;

Capp = zeros(size(Tguess));

if any(~tinyStep(:))
    flG = liquidFraction(Tguess(~tinyStep), Ts, Tl);
    flN = liquidFraction(Tn(~tinyStep),     Ts, Tl);
    Capp(~tinyStep) = rhocp + rhoL .* (flG - flN) ./ dT(~tinyStep);
end

if any(tinyStep(:))
    [~, dfldT] = liquidFraction(Tguess(tinyStep), Ts, Tl);
    Capp(tinyStep) = rhocp + rhoL .* dfldT;
end

end
