function matp = materialStruct(P)
%MATERIALSTRUCT  Pack the scalar material properties solveOneStepRadial
%needs (mushy-zone bounds + effective rho*cp / rho*L) from the parameter
%struct P. Kept as a tiny separate step so P itself stays the single
%source of truth (mirrors PCM_python/params.py's Params.material_dict()).

matp.Ts = P.Tm - P.dTmush/2;
matp.Tl = P.Tm + P.dTmush/2;
matp.kPCMs = P.kPCMs;
matp.kPCMl = P.kPCMl;
matp.rhocpPCM = P.rhoPCM * P.cpPCM;
matp.rhoLPCM  = P.rhoPCM * P.LPCM;
matp.kCu = P.kCu;
matp.rhocpCu = P.rhoCu * P.cpCu;
matp.Thotwater = P.Thotwater;

end
