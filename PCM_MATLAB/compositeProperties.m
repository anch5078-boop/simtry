function C = compositeProperties(P)
%COMPOSITEPROPERTIES  Effective (Cu-particle enhanced) PCM properties.
%
%   C = COMPOSITEPROPERTIES(P) computes the Maxwell effective
%   conductivity of the dispersed-copper-particle PCM composite,
%   separately for the solid and liquid base PCM conductivities, plus
%   the effective density and sensible/latent heat capacities. Because
%   the particle volume fraction phi is uniform and constant in the
%   Version-1 model, every field in C is a single scalar (evaluated
%   once, reused every time step).
%
%   Maxwell (Maxwell-Garnett) effective conductivity:
%
%       keff = km * (kCu + 2*km + 2*phi*(kCu-km)) ...
%                  / (kCu + 2*km -   phi*(kCu-km))
%
%   applied twice: once with km = kPCM_s (solid) and once with
%   km = kPCM_l (liquid).
%
%   Effective density (rule of mixtures):
%       rho_eff = (1-phi)*rhoPCM + phi*rhoCu
%
%   Effective sensible heat capacity:
%       (rho*cp)_eff = (1-phi)*rhoPCM*cpPCM + phi*rhoCu*cpCu
%
%   Effective latent heat capacity -- only the PCM fraction stores
%   latent heat, the copper particles do not melt:
%       (rho*L)_eff = (1-phi)*rhoPCM*LPCM

phi = P.phi;

C.ks_eff = maxwellKeff(P.kPCM_s, P.kCu, phi);
C.kl_eff = maxwellKeff(P.kPCM_l, P.kCu, phi);

C.rho_eff    = (1-phi)*P.rhoPCM + phi*P.rhoCu;
C.rhocp_eff  = (1-phi)*P.rhoPCM*P.cpPCM + phi*P.rhoCu*P.cpCu;
C.rhoL_eff   = (1-phi)*P.rhoPCM*P.LPCM;

% Heater (pure copper) apparent heat capacity -- no phase change.
C.rhocp_Cu = P.rhoCu * P.cpCu;
C.k_Cu     = P.kCu;

end

% ------------------------------------------------------------------------
function keff = maxwellKeff(km, kd, phi)
%MAXWELLKEFF  Maxwell (dilute-limit) effective conductivity of a
%continuous matrix km with spherical dispersed inclusions kd at volume
%fraction phi.
keff = km .* (kd + 2*km + 2*phi.*(kd-km)) ./ (kd + 2*km - phi.*(kd-km));
end
