function E = enhancementFactor(gap, b, dTchar, P)
%ENHANCEMENTFACTOR  E(gap): the engineering effective-conductivity
%enhancement factor used by radialEnthalpyModel.m for the melt layer.
%See that file's header comment for the full derivation/caveats.
%
%   E = ENHANCEMENTFACTOR(gap, b, dTchar, P)
%
%       gap    : current total melt-layer thickness [m]
%       b      : Navier slip length [m]
%       dTchar : characteristic driving temperature difference across
%                the melt layer (Th - Tm) [deg C]
%       P      : parameter struct (needs deltaC, g, betaLiq, nuLiq,
%                alphaLiq, Pr)
%
%   E = clip( max(slipTerm(gap), convTerm(gap)), 1, ECAP )
%
%   slipTerm(gap) = 1 + 3*(b/gap)*exp(-gap/deltaC)
%       -- decaying continuation of the flow-rate enhancement factor
%       derived in ccmSlipModel.m's Model A (1+3b/delta for b<<delta)
%       into the regime where the true dynamic-force-balance film no
%       longer applies (deltaC ~ a few hundred um sets where that decay
%       happens). This is a deliberate *approximation*: it borrows
%       Model A's flow-enhancement algebra as a proxy for near-wall
%       convective augmentation of heat transfer, which is a reasonable
%       qualitative stand-in but -- unlike Model A's force balance --
%       is not itself separately derived from first principles here.
%
%   convTerm(gap) = max(1, 0.386*(Pr/(0.861+Pr))^(1/4) * Ra_gap^(1/4))
%       -- standard Raithby-Hollands-type concentric-annulus natural
%       convection correlation, Ra_gap = g*betaLiq*dTchar*gap^3 /
%       (nuLiq*alphaLiq).
%
%   ECAP = 30 -- both terms are capped; near gap -> 0 the algebraic
%       forms are not meant to be extrapolated indefinitely (that
%       regime is where Model A's dynamic force balance actually
%       governs). 30x is in the range of reported CCM/slip enhancement
%       factors over plain conduction and is applied only as a
%       ceiling, rarely the binding value away from gap -> 0.
%
%   This is explicitly a *simplified engineering blend*, not a
%   spatially-resolved treatment of two distinct flow mechanisms.

ECAP = 30;

if gap <= 1e-9 || dTchar <= 0
    E = 1.0;
    return;
end

slipTerm = 1 + 3*(b/gap)*exp(-gap/P.deltaC);

Ra = P.g * P.betaLiq * dTchar * gap^3 / (P.nuLiq * P.alphaLiq);
convTerm = 1.0;
if Ra > 1e3
    Pr = P.Pr;
    convTerm = 0.386*(Pr/(0.861+Pr))^0.25 * Ra^0.25;
end

E = max([slipTerm, convTerm, 1.0]);
E = min(E, ECAP);

end
