function isOn = pulseState(t, P)
%PULSESTATE  True if the heater pulse is ON at time t [s]. Mirrors
%PCM_MATLAB/heaterPower.m's `mod(t, period) < ton` convention.
%
%   isOn = PULSESTATE(t, P)

if P.pulseOn
    period = P.tOn + P.tOff;
    isOn = mod(t, period) < P.tOn;
else
    isOn = true;
end

end
