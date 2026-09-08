"""
parameters.py -- central place for every physical, geometric and numerical
constant used by the Version-2 slip-enhanced close-contact-melting (sCCM)
model of a vertical-tube pulsed-heater PCM thermal battery.

Edit values here -- nothing else in the code should contain a "magic
number". Mirrors the single-params-file convention used in
``PCM_MATLAB/parameters.m`` for the Version-1 conduction-only model.

PLACEHOLDER DATA
-----------------
Geometry and PCM/heater thermophysical properties below are
representative placeholders (same paraffin-type PCM / copper heater
sleeve used as the Version-1 default, see ``PCM_MATLAB/parameters.m``)
picked to be self-consistent and to sit in a realistic range -- they are
NOT measured values for any specific device. Replace the block marked
"EDIT ME" with real datasheet / experimental numbers before drawing
quantitative conclusions. Two properties absent from the Version-1
conduction model are added here because the CCM/convection physics
needs them and they were not previously in the repo:

  * ``muL``   liquid dynamic viscosity  -- sets the CCM film thickness
              and the natural-convection Rayleigh number. This is the
              single most consequential unmeasured number in this
              model: CCM film thickness scales as muL**(1/4), so a 2x
              error in viscosity only shifts delta by ~19%, but
              viscosity itself can easily vary by 3-10x between PCM
              candidates and with temperature -- get a real value.
  * ``betaL`` liquid thermal expansion coefficient -- sets the
              buoyancy-driven natural-convection Rayleigh number.

Slip length ``b`` is set from the companion Nature paper this model is
grounded on (Li et al. 2026, "Pulse heating and slip enhance charging
of phase-change thermal batteries"), which reports a coating giving
b ~ 45-90 um via slip-enhanced close-contact melting (sCCM). Both ends
of that reported range are exposed here as named scenarios.
"""

from dataclasses import dataclass, field
from typing import Dict


@dataclass
class Params:
    # ---- Geometry (metres) --------------------------------------------
    # Vertical tube, radial melting: a cylindrical heater sleeve of
    # radius r_h runs up the axis of a tube of inner radius R; the PCM
    # occupies the annulus between them over the full column height
    # H_col, of which only the length H_heater (starting at z_base) is
    # actually heated. H_col >> H_heater so there is always a reservoir
    # of solid PCM above the heated section feeding the sCCM process
    # (Model A) -- a realistic "thermal battery" tube layout.
    r_h: float = 0.006       # heater sleeve outer radius            [m]
    R: float = 0.030         # container (tube) inner radius         [m]
    H_heater: float = 0.040  # heated length of the sleeve            [m]
    H_col: float = 0.150     # total PCM column height                [m]

    # ---- PCM material properties (paraffin-type, e.g. RT-58 class) ----
    # EDIT ME: replace with real datasheet values for your PCM.
    rho_s: float = 860.0      # solid density                        [kg/m^3]
    rho_l: float = 800.0      # liquid density                       [kg/m^3]
    cp: float = 2000.0        # specific heat (solid ~ liquid)       [J/(kg K)]
    k_s: float = 0.20         # solid conductivity                   [W/(m K)]
    k_l: float = 0.15         # liquid conductivity                  [W/(m K)]
    Lf: float = 180000.0      # latent heat of fusion                [J/kg]
    Ts: float = 55.0          # solidus temperature                  [deg C]
    Tl: float = 60.0          # liquidus temperature                 [deg C]
    mu_l: float = 0.010       # liquid dynamic viscosity  EDIT ME    [Pa s]
    beta_l: float = 8.0e-4    # liquid thermal expansion coefficient [1/K]

    # ---- Heater (pulsed wall boundary condition) -----------------------
    T_h: float = 90.0         # heater sleeve surface temp, ON        [deg C]
    T_init: float = 20.0      # initial PCM temperature, uniform      [deg C]
    pulse_on: bool = True     # False -> constant T_h (no pulsing)
    t_on: float = 30.0        # pulse ON duration                     [s]
    t_off: float = 30.0       # pulse OFF duration                    [s]

    # ---- Slip coating (grounding: Li et al. 2026 Nature sCCM paper) ---
    b_slip: float = 70.0e-6   # Navier slip length, default case      [m]
    delta0: float = 10.0e-6   # minimum achievable film thickness, a
                              # stand-in for surface roughness / other
                              # non-ideal-contact effects the smooth-
                              # wall lubrication force balance does not
                              # capture. Model A clamps its solved
                              # delta to this floor -- see README        [m]
    delta_c: float = 500.0e-6 # thin-film/slip-regime crossover scale
                              # used only by Model B's effective-
                              # conductivity submodel                  [m]

    # ---- Time stepping ---------------------------------------------------
    t_max: float = 7200.0     # safety cap on simulation time          [s]
    dt_A: float = 0.2         # Model A (0-D ODE) time step             [s]
    dt_B: float = 2.0         # Model B (1-D PDE) time step             [s]
    n_r: int = 240            # Model B radial control volumes          [-]

    g: float = 9.81           # gravitational acceleration            [m/s^2]

    # Named slip-length scenarios spanning the paper's reported range,
    # used by run_all.py for the comparison plots.
    slip_scenarios: Dict[str, float] = field(default_factory=lambda: {
        "no slip (b=0, baseline CCM)": 0.0,
        "slip b=45 um (paper, low end)": 45.0e-6,
        "slip b=90 um (paper, high end)": 90.0e-6,
    })

    # ---- Derived quantities ---------------------------------------------
    @property
    def Tm(self) -> float:
        """Nominal melting temperature (mushy-interval midpoint) used by
        the 0-D lumped Model A, which does not resolve a mushy zone."""
        return 0.5 * (self.Ts + self.Tl)

    @property
    def alpha_l(self) -> float:
        """Liquid thermal diffusivity [m^2/s]."""
        return self.k_l / (self.rho_l * self.cp)

    @property
    def nu_l(self) -> float:
        """Liquid kinematic viscosity [m^2/s]."""
        return self.mu_l / self.rho_l

    @property
    def Pr(self) -> float:
        """Liquid Prandtl number [-]."""
        return self.nu_l / self.alpha_l

    @property
    def A_cross(self) -> float:
        """Annular PCM cross-section area [m^2]."""
        import math
        return math.pi * (self.R ** 2 - self.r_h ** 2)

    def pulse_state(self, t: float) -> bool:
        """True if the heater pulse is ON at time t [s]. Mirrors
        PCM_MATLAB/heaterPower.m's ``mod(t, period) < ton`` convention."""
        if not self.pulse_on:
            return True
        period = self.t_on + self.t_off
        return (t % period) < self.t_on

    def wall_temp(self, t: float) -> float:
        """Heater sleeve surface temperature at time t: T_h during ON,
        and *undefined* during OFF -- Model B applies a zero-flux
        (adiabatic) condition during OFF instead of using this value;
        Model A likewise sets driving dT=0 during OFF. Provided for
        convenience/plotting only."""
        return self.T_h if self.pulse_state(t) else self.T_init
