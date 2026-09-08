function R = radialEnthalpyModel(P, b, saveEvery)
%RADIALENTHALPYMODEL  Model B: 1-D transient radial enthalpy-method
%solver for the PCM annulus surrounding a pulsed vertical heater
%sleeve, representative of a single axial cross-section through the
%heated section (fixed positions, no sinking -- contrast with Model A
%in ccmSlipModel.m).
%
%   R = RADIALENTHALPYMODEL(P, b) or R = RADIALENTHALPYMODEL(P, b, saveEvery)
%   runs the model for slip length b [m] and returns a struct R with
%   fields t, rC, rMelt, liquidFracVol, Eused, pulseOn, tSaved,
%   THistory (a cell array of temperature-profile snapshots), b,
%   milestones. saveEvery (default 30) controls how often a snapshot is
%   kept for THistory/tSaved.
%
%Gives what Model A cannot: an actual T(r,t) field and a melt-front
%rMelt(t) trajectory, including the visible effect of pulsing (melting
%stalls -- and the field relaxes -- during each OFF interval).
%
%NUMERICAL METHOD (mirrors PCM_MATLAB Version-1, adapted to 1-D
%cylindrical coordinates)
%--------------------------------------------------------------------
%Finite-volume, backward-Euler in time, harmonic-mean face
%conductivities, a sparse matrix assembled from (row,col,val) triplets
%exactly like PCM_MATLAB/solveOneTimeStep.m's assembleConductionMatrix
%(so `A \ rhs` -- no hand-rolled tridiagonal solver), a few Picard
%sub-iterations per step for the k(T)/Capp(T) nonlinearity, and the
%same *secant* (energy-conserving) apparent heat capacity used in
%apparentCapacity.m / PCM_MATLAB/apparentCapacity.m. The domain is
%initialized at P.Tinit (well below Ts), not at Tm, for the same reason
%documented in PCM_MATLAB/README.md (initializing exactly at Tm falsely
%puts every cell at fl=0.5 before it has absorbed any real heat).
%
%Inner boundary at r=rH: Dirichlet T=Th while the pulse is ON, and
%adiabatic (no term added -- true zero flux, not just "no change")
%while OFF. Outer boundary at r=R: always adiabatic (insulated
%container wall), matching the Version-1 model's outer-boundary
%convention.
%
%EFFECTIVE-CONDUCTIVITY SUBMODEL FOR THE MELT LAYER
%------------------------------------------------------
%Resolving the sub-mm slip-lubricated film and the natural-convection
%circulation cell directly on this 1-D grid is out of scope (that is a
%multi-dimensional flow problem); instead, once PCM at a location has
%melted, its *conductivity* is scaled up by an engineering enhancement
%factor E(gap) (see enhancementFactor.m) intended to reproduce the
%right order of magnitude and the right qualitative regime crossover,
%where gap = the current total melt-layer thickness, not each cell's
%own depth -- this mirrors how natural-convection "effective
%conductivity" correlations (e.g. Raithby & Hollands 1975 for
%concentric-cylinder annuli) are normally used: one k_eff for the whole
%gap, replacing the conduction profile, not a locally-varying value
%within it. This is explicitly a *simplified engineering blend*, not a
%spatially-resolved treatment of two distinct flow mechanisms.
%
%A NOTE ON TIME LABELING
%------------------------
%Each step advances the field from t_old = (step-1)*dt to
%t_new = t_old + dt using a boundary condition decided by
%pulseState(t_old, P). Every quantity recorded for this step (T,
%rMelt, liquidFracVol, ...) describes the field at t_new and is stored
%against that time -- not t_old -- since that is the time the recorded
%state actually corresponds to.

if nargin < 3
    saveEvery = 30;
end

N  = P.nR;
dr = (P.R - P.rH) / N;
rFace = P.rH + (0:N)'*dr;
rC    = P.rH + ((0:N-1)' + 0.5)*dr;
V     = pi*(rFace(2:end).^2 - rFace(1:end-1).^2);   % per unit height

T  = P.Tinit * ones(N,1);
dt = P.dtB;
nSteps = floor(P.tMax/dt) + 1;

tList       = zeros(nSteps,1);
rMeltList   = zeros(nSteps,1);
liqFracList = zeros(nSteps,1);
Elist       = zeros(nSteps,1);
pulseList   = false(nSteps,1);

tSavedList = [];
THistList  = {};

dTchar = P.Th - P.Tm;
picardIters = 4;
picardTol   = 1e-3;

p = (1:N-1)';   % left cell of each interior face (used every Picard iterate)
q = (2:N)';     % right cell of each interior face

nUsed = nSteps;
for step = 1:nSteps
    tOld = (step-1)*dt;   % computed fresh each step (not accumulated) to
                          % avoid float drift landing on the wrong side
                          % of a pulse edge
    tNew = tOld + dt;
    isOn = pulseState(tOld, P);   % BC applied while advancing tOld -> tNew

    % Current melt-front / gap estimate (from the previous step's
    % field) drives this step's enhancement factor; re-evaluated once
    % more after the implicit solve below (one lag is fine -- gap
    % changes slowly compared to dtB).
    flPrev = liquidFraction(T, P.Ts, P.Tl);
    gap = dr * sum(flPrev);
    if isOn
        E = enhancementFactor(gap, b, dTchar, P);
    else
        E = 1.0;
    end

    Tn     = T;
    Tguess = T;
    for pit = 1:picardIters
        fl = liquidFraction(Tguess, P.Ts, P.Tl);
        k  = (1-fl)*P.kSol + fl*(P.kLiq*E);
        Capp = apparentCapacity(Tguess, Tn, P.Ts, P.Tl, P.rhoSol*P.cp, P.rhoSol*P.Lf);

        kFace  = 2*k(1:end-1).*k(2:end) ./ (k(1:end-1) + k(2:end));   % interior faces, N-1 of them
        GcFace = kFace .* 2*pi.*rFace(2:end-1) / dr;                   % interior conductances

        capTerm = Capp .* V / dt;

        % Sparse assembly, same (row,col,val) triplet pattern as
        % PCM_MATLAB/solveOneTimeStep.m's assembleConductionMatrix:
        % each interior face contributes -Gc off-diagonal and +Gc onto
        % both neighbouring diagonals (duplicate (row,col) triplets are
        % summed by `sparse`, which is how the capacity-diagonal
        % triplets below combine with the face contributions).
        rows = [p; q; p; q];
        cols = [q; p; p; q];
        vals = [-GcFace; -GcFace; GcFace; GcFace];

        rows = [rows; (1:N)'];
        cols = [cols; (1:N)'];
        vals = [vals; capTerm];

        rhs = capTerm .* Tn;

        % inner wall (r=rH): Dirichlet Th while ON, adiabatic while OFF
        if isOn
            GcWall = k(1) * 2*pi*P.rH / (0.5*dr);
            rows = [rows; 1];
            cols = [cols; 1];
            vals = [vals; GcWall];
            rhs(1) = rhs(1) + GcWall*P.Th;
        end
        % outer wall (r=R): always adiabatic -- no term added

        A = sparse(rows, cols, vals, N, N);
        Tnew = A \ rhs;

        change = max(abs(Tnew - Tguess));
        scale  = max(1, max(abs(Tguess)));
        Tguess = Tnew;
        if change < picardTol*scale
            break;
        end
    end

    T = Tguess;
    flNow = liquidFraction(T, P.Ts, P.Tl);

    rMelt = frontPosition(rC, T, P.Tm);
    liqFracVol = sum(flNow .* V) / sum(V);

    tList(step)       = tNew;
    rMeltList(step)   = rMelt;
    liqFracList(step) = liqFracVol;
    Elist(step)       = E;
    pulseList(step)   = isOn;

    if mod(step-1, saveEvery) == 0
        THistList{end+1} = T;   %#ok<AGROW>
        tSavedList(end+1) = tNew;   %#ok<AGROW>
    end

    if liqFracVol > 0.999
        nUsed = step;
        break;
    end
end

R.t             = tList(1:nUsed);
R.rC            = rC;
R.rMelt         = rMeltList(1:nUsed);
R.liquidFracVol = liqFracList(1:nUsed);
R.Eused         = Elist(1:nUsed);
R.pulseOn       = pulseList(1:nUsed);
R.tSaved        = tSavedList(:);
R.THistory      = THistList;
R.b             = b;
R.milestones    = milestonesFromSeries(R.t, R.liquidFracVol);

end
