function Capp = apparentCapacity(Tguess, Tn, Ts, Tl, rhocp, rhoL)
%APPARENTCAPACITY  Energy-conserving ("secant") apparent heat capacity
%for the phase-change term, used to build both the diagonal and the
%right-hand-side of the implicit conduction system with a *consistent*
%capacity value so the discrete equation exactly conserves the
%sensible + latent enthalpy change between the previous time level Tn
%and the current Picard iterate Tguess:
%
%   rhocp*(Tguess-Tn) + rhoL*(fl(Tguess)-fl(Tn)) = Capp*(Tguess-Tn)
%
%   => Capp = rhocp + rhoL * (fl(Tguess)-fl(Tn)) / (Tguess-Tn)
%
%This "chord slope" (rather than the local tangent dfl/dT) is what
%makes the scheme robust when a single time step's temperature change
%is comparable to or larger than the width of the smoothed mushy
%interval (Tl-Ts): the plain tangent-derivative apparent-heat-capacity
%method can then "jump across" the latent-heat spike within one Picard
%iterate and silently under-count the latent heat, making the melt
%front advance too fast (this is exactly the failure mode caught by
%validateStefan.m). The secant form captures the correct total latent
%heat absorbed no matter how large the step is.
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
