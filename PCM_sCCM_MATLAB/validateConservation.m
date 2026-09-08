%VALIDATECONSERVATION  Energy-conservation check for Model B
%(radialEnthalpyModel.m), in the same spirit as
%PCM_MATLAB/validateStefan.m's validation of the Version-1 solver: this
%codebase should hold itself to the same standard.
%
%What it checks: with both boundaries adiabatic during a pulse OFF
%interval (no source, no sink anywhere in the domain), the total
%sensible+latent enthalpy
%
%   H(t) = sum_i [ rhoSol*cp*T_i(t) + rhoSol*Lf*fl_i(T_i(t)) ] * V_i
%
%must stay constant between two time points that are both within the
%same OFF interval, to within time-discretization/Picard-convergence
%error. This directly exercises the exact defect this script (and its
%Python counterpart, validate_conservation.py, in the companion
%PCM_sCCM/ folder) was written to catch during development: an earlier
%version of the underlying solve had a sign error in how the
%conductance terms were assembled, which still produced a
%converged-looking Picard iteration but silently violated conservation
%-- H measurably (not just to floating point) drained away during every
%OFF interval, with no error or warning anywhere else. That failure
%mode would not have been caught by eyeballing T(r,t) plots alone. This
%MATLAB port avoids that whole bug class by assembling the conduction
%matrix with the same (row,col,val) sparse-triplet pattern as
%PCM_MATLAB/solveOneTimeStep.m (verified against a dense reference
%matrix during development) rather than a hand-rolled tridiagonal
%solver, but this check is kept as a standing regression guard in case
%the assembly is ever changed.
%
%Usage:
%   validateConservation
%Prints PASS/FAIL, mirroring validateStefan.m's report style.

REL_TOL = 0.02;   % max allowed |dH| during an OFF step, relative to the
                   % typical |dH| magnitude seen during an ON step in
                   % the same run (a "how much does this look like a
                   % real heat-input step" scale, not an absolute number)

P = parameters();
P.tMax = 200;   % short run (a few pulse cycles) -- this check only
                % needs the early transient, not a full charge
b = 45.0e-6;

R = radialEnthalpyModel(P, b, 1);   % saveEvery=1: keep every step's snapshot

dr = R.rC(2) - R.rC(1);
rFace = P.rH + (0:numel(R.rC))'*dr;
V = pi*(rFace(2:end).^2 - rFace(1:end-1).^2);

nSaved = numel(R.tSaved);
H = zeros(nSaved,1);
for i = 1:nSaved
    T = R.THistory{i};
    fl = liquidFraction(T, P.Ts, P.Tl);
    H(i) = sum((P.rhoSol*P.cp*T + P.rhoSol*P.Lf*fl) .* V);
end

% saveEvery=1 here, so tSaved/THistory line up 1:1 with the per-step
% pulseOn record. dH(i) = H(i+1)-H(i) is produced by the solve that
% generated THistory{i+1}, whose boundary condition is pulseOn(i+1) --
% so index with pulseOn(2:end), not pulseOn(1:end-1), either of which
% would mislabel every step spanning an ON/OFF transition.
isOn = R.pulseOn(2:nSaved);
dH = diff(H);

if ~any(isOn) || ~any(~isOn)
    error('validateConservation:noTransition', ...
        'FAIL: run did not include both an ON and an OFF step -- check parameters tOn/tOff/tMax.');
end

scale = mean(abs(dH(isOn)));
offDrift = max(abs(dH(~isOn)));
relDrift = offDrift / scale;

fprintf('Typical |dH| per step while ON  : %10.4f  (energy actually entering)\n', scale);
fprintf('Worst   |dH| per step while OFF : %10.4f  (should be ~0 -- adiabatic)\n', offDrift);
fprintf('Ratio (worst OFF drift / typical ON step): %.4f%%\n', 100*relDrift);

if relDrift < REL_TOL
    fprintf('PASS: OFF-interval energy drift is %.4f%% of a typical ON-step input, well under the %.0f%% tolerance -- Model B conserves energy.\n', ...
        100*relDrift, 100*REL_TOL);
else
    error('validateConservation:driftTooLarge', ...
        'FAIL: OFF-interval energy drift is %.4f%% of a typical ON-step input, exceeding the %.0f%% tolerance.', ...
        100*relDrift, 100*REL_TOL);
end
