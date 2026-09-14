function [Tnew, fl] = solveOneStepRadial(Tn, G, matp, Qtotal, dt, picardIters, picardTol)
%SOLVEONESTEPRADIAL  Advance the radial temperature field by one implicit
%(backward-Euler) time step on the grid G, solving
%
%   Capp(T) * V/dt * (T - Tn) = sum_faces G_face*(T_nb - T) + Q_cell
%
%with a Dirichlet inner boundary (fixed at matp.Thotwater, the pipe
%surface) and an adiabatic outer boundary (the insulated container
%wall, reproduced simply by never adding a flux term there). k(T) and
%Capp(T) are nonlinear through the liquid fraction, so a few Picard
%sub-iterations are used, exactly as in PCM_MATLAB/solveOneTimeStep.m.
%
%   [Tnew, fl] = SOLVEONESTEPRADIAL(Tn, G, matp, Qtotal, dt, picardIters, picardTol)
%
%       Tn      : N x 1 temperature at the start of the step   [deg C]
%       G       : geometry struct from createGeometryRadial.m / createGeometryNoSleeve.m
%       matp    : material struct from materialStruct.m
%       Qtotal  : instantaneous sleeve power this step          [W/m of pipe length]
%       dt      : time step                                     [s]
%
%       Tnew, fl : N x 1 temperature / liquid fraction at the end of the step

N = G.N;
V = G.V;
dr = G.dr;
isCu  = G.isCu;
isPcm = G.isPcm;

VCu = sum(V(isCu));
qCell = zeros(N, 1);
if VCu > 0 && Qtotal ~= 0
    qCell(isCu) = Qtotal * V(isCu) / VCu;   % distribute by volume
end

Tguess = Tn;
idx = (1:N)';

for it = 1:picardIters
    fl = liquidFractionRadial(Tguess, matp.Ts, matp.Tl);

    k = zeros(N, 1);
    k(isPcm) = (1 - fl(isPcm)).*matp.kPCMs + fl(isPcm).*matp.kPCMl;
    k(isCu)  = matp.kCu;

    Capp = zeros(N, 1);
    Capp(isPcm) = apparentCapacityRadial(Tguess(isPcm), Tn(isPcm), ...
                                          matp.Ts, matp.Tl, matp.rhocpPCM, matp.rhoLPCM);
    Capp(isCu)  = matp.rhocpCu;

    % Harmonic-mean interior face conductance (uniform dr -> cell
    % centres equidistant from the shared face).
    kFace = 2*k(1:end-1).*k(2:end) ./ (k(1:end-1) + k(2:end));
    Gint  = kFace .* G.Aface(2:end-1) / dr;         % N-1 x 1

    % Inner Dirichlet boundary conductance (pipe surface, dr/2 from the
    % first cell centre).
    GbcIn = 2*k(1)*G.Aface(1) / dr;

    diagMain = Capp .* V / dt;
    diagMain(1:end-1) = diagMain(1:end-1) + Gint;
    diagMain(2:end)   = diagMain(2:end)   + Gint;
    diagMain(1)       = diagMain(1) + GbcIn;

    offDiag = -Gint;

    rows = [idx(1:end-1); idx(2:end); idx];
    cols = [idx(2:end); idx(1:end-1); idx];
    vals = [offDiag; offDiag; diagMain];
    A = sparse(rows, cols, vals, N, N);

    b = Capp .* V / dt .* Tn + qCell;
    b(1) = b(1) + GbcIn * matp.Thotwater;

    Tnew = A \ b;

    change = max(abs(Tnew - Tguess));
    scale  = max(1, max(abs(Tguess)));
    Tguess = Tnew;

    if change < picardTol * scale
        break;
    end
end

fl = liquidFractionRadial(Tnew, matp.Ts, matp.Tl);

end
