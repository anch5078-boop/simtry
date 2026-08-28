function P = parameters()
%PARAMETERS  Central place for every physical, geometric and numerical
%constant used by the Version-1 (conduction + phase change) PCM solver.
%
%   P = PARAMETERS() returns a struct with all inputs the rest of the
%   code needs. Edit values here -- nothing else in the code should
%   contain a "magic number".
%
%   Default material values are representative of a paraffin-type PCM
%   (e.g. RT-58 class wax) heated by two copper cylindrical sleeves.
%   Replace with datasheet values for your actual PCM before publishing
%   results.

%% ---- Geometry (metres) --------------------------------------------
P.R        = 0.030;   % outer PCM cylinder radius            [m]
P.H        = 0.060;   % cylinder height                      [m]

P.r1i      = 0.006;   % heater 1 (inner) inner radius         [m]
P.th1      = 0.002;   % heater 1 wall thickness               [m]
P.r1o      = P.r1i + P.th1;

P.r2i      = 0.018;   % heater 2 (outer) inner radius         [m]
P.th2      = 0.002;   % heater 2 wall thickness               [m]
P.r2o      = P.r2i + P.th2;

P.zbase    = 0.010;   % heater vertical start                 [m]
P.Hheater  = 0.040;   % heater vertical extent                [m]

%% ---- Grid ------------------------------------------------------------
% Cubic finite-volume cells of size h. Domain is a Cartesian box
% [-R,R] x [-R,R] x [0,H] with only cells satisfying x^2+y^2<=R^2
% (and 0<=z<=H) marked "active". Everything outside the cylinder is
% simply not simulated (equivalent to a perfectly insulated outer
% cylindrical wall).
P.h        = 0.0015;  % cell size                             [m]

%% ---- PCM material properties ------------------------------------
P.rhoPCM   = 800;      % PCM density                          [kg/m^3]
P.cpPCM    = 2000;     % PCM specific heat                    [J/(kg K)]
P.kPCM_s   = 0.20;     % PCM conductivity, solid phase        [W/(m K)]
P.kPCM_l   = 0.15;     % PCM conductivity, liquid phase       [W/(m K)]
P.LPCM     = 180000;   % PCM latent heat of fusion            [J/kg]
P.Ts       = 55;       % solidus temperature                  [deg C]
P.Tl       = 60;       % liquidus temperature                 [deg C]

%% ---- Copper (heater sleeves + dispersed particles) ----------------
P.rhoCu    = 8960;     % copper density                       [kg/m^3]
P.cpCu     = 385;      % copper specific heat                 [J/(kg K)]
P.kCu      = 400;      % copper conductivity                  [W/(m K)]

%% ---- Dispersed copper particle enhancement -----------------------
P.phi      = 0.00;     % Cu particle volume fraction in PCM   [-]
                        % (0 for pure-PCM cases C0/C1)

%% ---- Heater power / pulsing ---------------------------------------
P.Ppeak      = 8.0;    % total peak heater power              [W]
P.pulseOn    = false;  % false -> constant power = Ppeak
P.ton        = 30;     % pulse ON duration                    [s]
P.toff       = 30;     % pulse OFF duration                   [s]

%% ---- Initial / reference conditions --------------------------------
P.Tinit    = 20;       % initial temperature, uniform         [deg C]

%% ---- Time stepping --------------------------------------------------
P.dt           = 0.5;      % time step                        [s]
P.tMax         = 3600;     % safety cap on simulation time     [s]
P.picardIters  = 2;        % nonlinear (Picard) sub-iterations per step
P.picardTol    = 1e-3;     % relative T tolerance to stop Picard loop early

%% ---- Output control --------------------------------------------------
P.saveEvery    = 20;    % store a snapshot every N time steps

end
