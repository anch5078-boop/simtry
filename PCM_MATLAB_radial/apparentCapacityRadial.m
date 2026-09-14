function Capp = apparentCapacityRadial(Tguess, Tn, Ts, Tl, rhocp, rhoL)
%APPARENTCAPACITYRADIAL  Energy-conserving ("secant") apparent heat
%capacity -- identical in form to PCM_MATLAB/apparentCapacity.m (see
%that file's header comment for the full derivation and for why the
%naive tangent-derivative form under-counts latent heat when a step's
%dT is comparable to the mushy-zone width):
%
%   Capp = rhocp + rhoL * (fl(Tguess)-fl(Tn)) / (Tguess-Tn)
%
%falling back to the tangent derivative dfl/dT where Tguess and Tn
%coincide to machine precision.

dT = Tguess - Tn;
tinyStep = abs(dT) < 1e-8;

Capp = zeros(size(Tguess));

if any(~tinyStep(:))
    flG = liquidFractionRadial(Tguess(~tinyStep), Ts, Tl);
    flN = liquidFractionRadial(Tn(~tinyStep),     Ts, Tl);
    Capp(~tinyStep) = rhocp + rhoL .* (flG - flN) ./ dT(~tinyStep);
end

if any(tinyStep(:))
    dfldT = liquidFractionDerivRadial(Tguess(tinyStep), Ts, Tl);
    Capp(tinyStep) = rhocp + rhoL .* dfldT;
end

end
