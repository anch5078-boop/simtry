function [Qv, P1, P2, pulseOn] = heaterPower(t, G, P)
%HEATERPOWER  Volumetric heat generation Q''' [W/m^3] in the two
%cylindrical heater sleeves at time t, including optional square-wave
%pulsing. The two heaters always share the instantaneous total power
%equally (P1 = P2 = Ptot/2).
%
%   [Qv, P1, P2, pulseOn] = HEATERPOWER(t, G, P)
%
%       Qv      : Na x 1 volumetric heat generation, nonzero only in
%                 heater-1 / heater-2 cells (0 elsewhere, incl. PCM)
%       P1, P2  : instantaneous power [W] delivered by each heater
%       pulseOn : logical, true if the pulse is currently ON
%                 (always true when P.pulseOn == false, i.e. constant
%                 power mode)

if P.pulseOn
    period  = P.ton + P.toff;
    pulseOn = mod(t, period) < P.ton;
else
    pulseOn = true;
end

if pulseOn
    Ptotal = P.Ppeak;
else
    Ptotal = 0;
end

P1 = Ptotal / 2;
P2 = Ptotal / 2;

Qv = zeros(G.Na, 1);
Qv(G.heater1MaskA) = P1 / G.Vh1;
Qv(G.heater2MaskA) = P2 / G.Vh2;

end
