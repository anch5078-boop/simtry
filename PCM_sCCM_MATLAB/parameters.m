function P = parameters()
%PARAMETERS  Central place for every physical, geometric and numerical
%constant used by the Version-2 slip-enhanced close-contact-melting
%(sCCM) model of a vertical-tube pulsed-heater PCM thermal battery.
%
%   P = PARAMETERS() returns a struct with all inputs the rest of the
%   code needs. Edit values here -- nothing else in the code should
%   contain a "magic number". Mirrors the single-params-file
%   convention used by PCM_MATLAB/parameters.m for the Version-1
%   conduction-only model.
%
%PLACEHOLDER DATA
%-----------------
%Geometry and PCM/heater thermophysical properties below are
%representative placeholders (same paraffin-type PCM / copper heater
%sleeve used as the Version-1 default) picked to be self-consistent and
%to sit in a realistic range -- they are NOT measured values for any
%specific device. Replace before drawing quantitative conclusions.
%
%Two properties are new versus the Version-1 conduction-only model
%because the CCM/convection physics needs them:
%
%   muLiq   liquid dynamic viscosity -- sets the CCM film thickness
%           and the natural-convection Rayleigh number. This is the
%           single most consequential unmeasured number in this model:
%           CCM film thickness scales as muLiq^(1/4), so a 2x error in
%           viscosity only shifts delta by ~19%, but viscosity itself
%           can easily vary by 3-10x between PCM candidates and with
%           temperature -- get a real value.
%   betaLiq liquid thermal expansion coefficient -- sets the
%           buoyancy-driven natural-convection Rayleigh number.
%
%Slip length bSlip is grounded on Li et al. 2026 (Nature), "Pulse
%heating and slip enhance charging of phase-change thermal batteries",
%which reports a coating giving a Navier slip length b ~ 45-90 um via
%slip-enhanced close-contact melting (sCCM). runAll.m sweeps that whole
%reported range; bSlip here is just the single-case default.

%% ---- Geometry (metres) ----------------------------------------------
% Vertical tube, radial melting: a cylindrical heater sleeve of radius
% rH runs up the axis of a tube of inner radius R; the PCM occupies the
% annulus between them over the full column height Hcol, of which only
% the length Hheater (at the bottom of the column) is actually heated.
% Hcol >> Hheater so there is always a reservoir of solid PCM above the
% heated section feeding the sCCM process (Model A) -- a realistic
% "thermal battery" tube layout.
P.rH        = 0.006;    % heater sleeve outer radius             [m]
P.R         = 0.030;    % container (tube) inner radius          [m]
P.Hheater   = 0.040;    % heated length of the sleeve             [m]
P.Hcol      = 0.150;    % total PCM column height                 [m]

%% ---- PCM material properties (paraffin-type, e.g. RT-58 class) ------
% EDIT ME: replace with real datasheet values for your PCM.
P.rhoSol    = 860;      % solid density                          [kg/m^3]
P.rhoLiq    = 800;      % liquid density                         [kg/m^3]
P.cp        = 2000;     % specific heat (solid ~ liquid)         [J/(kg K)]
P.kSol      = 0.20;     % solid conductivity                     [W/(m K)]
P.kLiq      = 0.15;     % liquid conductivity                    [W/(m K)]
P.Lf        = 180000;   % latent heat of fusion                  [J/kg]
P.Ts        = 55;       % solidus temperature                    [deg C]
P.Tl        = 60;       % liquidus temperature                   [deg C]
P.muLiq     = 0.010;    % liquid dynamic viscosity  EDIT ME      [Pa s]
P.betaLiq   = 8.0e-4;   % liquid thermal expansion coefficient   [1/K]

%% ---- Heater (pulsed wall boundary condition) -------------------------
P.Th        = 90;       % heater sleeve surface temp, ON          [deg C]
P.Tinit     = 20;       % initial PCM temperature, uniform        [deg C]
P.pulseOn   = true;     % false -> constant Th (no pulsing)
P.tOn       = 30;       % pulse ON duration                       [s]
P.tOff      = 30;       % pulse OFF duration                      [s]

%% ---- Slip coating (grounding: Li et al. 2026 Nature sCCM paper) ------
P.bSlip     = 70.0e-6;  % Navier slip length, single-case default [m]
P.delta0    = 10.0e-6;  % minimum achievable film thickness, a
                         % stand-in for surface roughness / other
                         % non-ideal-contact effects the smooth-wall
                         % lubrication force balance does not capture.
                         % Model A clamps its solved delta to this
                         % floor -- see ccmSlipModel.m README notes.  [m]
P.deltaC    = 500.0e-6; % thin-film/slip-regime crossover scale used
                         % only by Model B's effective-conductivity
                         % submodel                                  [m]

%% ---- Time stepping -----------------------------------------------------
P.tMax      = 7200;     % safety cap on simulation time            [s]
P.dtA       = 1.0;      % Model A (0-D ODE) time step               [s]
                         % (0.2s matches the companion Python model
                         % bit-for-bit closer, but costs ~5x the
                         % runtime for a <0.01% milestone difference in
                         % Octave's interpreted loop -- see README)
P.dtB       = 2.0;      % Model B (1-D PDE) time step               [s]
P.nR        = 240;      % Model B radial control volumes            [-]

P.g         = 9.81;     % gravitational acceleration              [m/s^2]

%% ---- Derived quantities ------------------------------------------------
P.Tm       = 0.5*(P.Ts + P.Tl);              % nominal melting temperature
                                              % (mushy-interval midpoint)
                                              % used by lumped Model A
P.alphaLiq = P.kLiq / (P.rhoLiq*P.cp);       % liquid thermal diffusivity [m^2/s]
P.nuLiq    = P.muLiq / P.rhoLiq;             % liquid kinematic viscosity [m^2/s]
P.Pr       = P.nuLiq / P.alphaLiq;           % liquid Prandtl number [-]
P.Across   = pi*(P.R^2 - P.rH^2);            % annular PCM cross-section [m^2]

end
