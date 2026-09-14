function P = parametersRadial()
%PARAMETERSRADIAL  Every physical, geometric and control constant for the
%"hot pipe / PCM / pulsed copper sleeve / PCM / insulated wall" case, as
%drawn (sleeve centered in the pipe->wall gap). This is the MATLAB/Octave
%counterpart of PCM_python/params.py -- see that file's comments for the
%reasoning behind each assumed value; nothing here should need to change
%except by editing this file.

%% ---- Geometry (metres) -- rPipe is given, the rest is ASSUMED --------
P.rPipe   = 0.0127;    % 1" pipe OD / PCM-facing radius, given         [m]
P.rWall   = 0.045;     % container inner radius (ASSUMED)              [m]
P.tSleeve = 0.0015;    % copper sleeve radial thickness (ASSUMED)      [m]
P.centerFrac = 0.5;    % sleeve mid-radius as a fraction of the pipe->
                       % wall gap; 0.5 = "centered in the gap" as drawn

%% ---- Grid -------------------------------------------------------------
P.nCells = 240;

%% ---- PCM properties (paraffin/RT42-class wax, Tm ~ 40 C) -- ASSUMED,
%% replace with your PCM's real datasheet values.
P.rhoPCM = 830;        % kg/m^3 (single value, both phases)
P.cpPCM  = 2100;       % J/(kg K)
P.kPCMs  = 0.20;       % W/(m K), solid
P.kPCMl  = 0.15;       % W/(m K), liquid
P.LPCM   = 170000;     % J/kg, latent heat of fusion
P.Tm     = 40;         % deg C, nominal melt point -- given
P.dTmush = 4;          % deg C, smoothed mushy-zone width (numerical
                       % regularization only)

%% ---- Copper (sleeve) ----------------------------------------------------
P.rhoCu = 8960;
P.cpCu  = 385;
P.kCu   = 400;

%% ---- Boundary / operating conditions -- given --------------------------
P.Thotwater = 120;     % deg C, Dirichlet at the pipe surface (pipe wall
                       % resistance neglected)
P.Tinit     = 20;      % deg C, uniform initial temperature (ASSUMED)

%% ---- Sleeve pulse-heating control ---------------------------------------
% "Pulse heated at lowest power to reach 40 C": bisect for the smallest
% sleeve peak power that reaches Tm within tRampTarget seconds of a cold
% start, then hold there with a bang-bang thermostat (see
% runFullMeltRadial.m).
P.tRampTarget = 300;   % s (ASSUMED design target)

%% ---- Time stepping -------------------------------------------------------
P.dt = 2;              % s
P.tMax = 6*3600;       % s, safety cap
P.picardIters = 4;
P.picardTol   = 1e-3;

end
