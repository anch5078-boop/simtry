function G = createGeometryNoSleeve(P)
%CREATEGEOMETRYNOSLEEVE  Same radial grid as createGeometryRadial.m but
%with no copper sleeve anywhere (pure PCM from pipe to wall) -- the
%"pipe conduction only" reference case.

N = P.nCells;
rFace = linspace(P.rPipe, P.rWall, N+1)';
rc    = 0.5*(rFace(1:end-1) + rFace(2:end));
dr    = rFace(2) - rFace(1);
V     = pi*(rFace(2:end).^2 - rFace(1:end-1).^2);
Aface = 2*pi*rFace;

G.N = N;
G.rFace = rFace;
G.rc = rc;
G.dr = dr;
G.V = V;
G.Aface = Aface;
G.isCu  = false(N,1);
G.isPcm = true(N,1);
G.rSleeveIn = NaN;
G.rSleeveOut = NaN;

end
