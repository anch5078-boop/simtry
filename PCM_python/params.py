"""
params.py -- every physical, geometric and control constant for the
"hot pipe / PCM / pulsed copper sleeve / PCM / insulated wall" case.

Values not given in the task prompt (container radius, sleeve
thickness, PCM material dataset, ramp-up target, grid resolution) are
engineering defaults, clearly flagged below -- edit them here to match
an actual PCM datasheet / hardware design; nothing else in the code
should need to change.
"""

from dataclasses import dataclass, field


@dataclass
class Params:
    # ---- Geometry (metres) -- from the prompt -----------------------
    r_pipe: float = 0.0127          # 1" pipe OD / PCM-facing radius, given [m]

    # ASSUMED (not given): overall container size and sleeve thickness.
    # "Wider gap: pipe to wall" in the sketch -> pick a lab/prototype-scale
    # annulus a few pipe-radii wide.
    r_wall: float = 0.045           # container inner radius            [m]  (ASSUMED)
    t_sleeve: float = 0.0015        # copper sleeve radial thickness     [m]  (ASSUMED, thin)
    sleeve_center_frac: float = 0.5  # sleeve mid-radius as a fraction of
                                     # the pipe->wall gap; 0.5 = "centered
                                     # in the gap" as drawn.              (ASSUMED default)

    # ---- Grid --------------------------------------------------------
    n_cells: int = 240

    # ---- PCM properties (paraffin-wax class material, Tm ~ 40 C, e.g.
    # RT42-like commercial wax) -- ASSUMED representative values; replace
    # with the real datasheet before publishing results.
    rho_pcm: float = 830.0          # kg/m^3 (single value for both phases,
                                     # matching the modelling convention
                                     # already used in PCM_MATLAB/)
    cp_pcm: float = 2100.0          # J/(kg K)
    k_pcm_s: float = 0.20           # W/(m K), solid
    k_pcm_l: float = 0.15           # W/(m K), liquid
    L_pcm: float = 170000.0         # J/kg, latent heat of fusion
    Tm: float = 40.0                # deg C, nominal melt point -- given
    dT_mush: float = 4.0            # deg C, smoothed mushy-zone width (numerical
                                     # regularization only; Ts=Tm-dT/2, Tl=Tm+dT/2)

    # ---- Copper (sleeve) ----------------------------------------------
    rho_cu: float = 8960.0
    cp_cu: float = 385.0
    k_cu: float = 400.0

    # ---- Boundary / operating conditions -- given -----------------------
    T_hotwater: float = 120.0       # deg C, Dirichlet at pipe surface (pipe
                                     # wall resistance neglected -- see
                                     # model.py docstring)
    T_init: float = 20.0            # deg C, uniform initial temperature (ASSUMED)

    # ---- Sleeve pulse-heating control -----------------------------------
    # "Pulse heated at lowest power to reach 40 C": the sleeve is driven at
    # a fixed peak power P_peak whenever it is below Tm (full-on ramp),
    # and switched off once it reaches Tm (bang-bang thermostat holding it
    # at the melting point thereafter -- the standard idealisation of a
    # feedback-controlled resistive/induction heater, and a good
    # approximation for a thin, near-isothermal copper sleeve). P_peak
    # itself is *found*, not assumed: see sweep.find_min_power -- it is
    # the smallest peak power for which the sleeve reaches Tm within
    # `t_ramp_target` seconds of start-up.
    t_ramp_target: float = 300.0    # s, target sleeve ramp-up time (ASSUMED
                                     # design target -- "pulse heated
                                     # quickly to the melt point", not a
                                     # slow multi-hour creep)

    # ---- Time stepping ---------------------------------------------------
    dt: float = 2.0                 # s
    t_max: float = 6 * 3600.0       # s, safety cap
    picard_iters: int = 4
    picard_tol: float = 1e-3

    def mush_bounds(self):
        return self.Tm - self.dT_mush / 2, self.Tm + self.dT_mush / 2

    def material_dict(self):
        Ts, Tl = self.mush_bounds()
        return dict(
            Ts=Ts, Tl=Tl,
            k_pcm_s=self.k_pcm_s, k_pcm_l=self.k_pcm_l,
            rhocp_pcm=self.rho_pcm * self.cp_pcm,
            rhoL_pcm=self.rho_pcm * self.L_pcm,
            k_cu=self.k_cu, rhocp_cu=self.rho_cu * self.cp_cu,
            T_hotwater=self.T_hotwater,
        )

    def sleeve_radii(self, center_frac=None):
        """Inner/outer radius of the copper sleeve for a given fractional
        position across the pipe->wall gap (0 = touching pipe, 1 =
        touching wall). Defaults to sleeve_center_frac ("centered in the
        gap" per the diagram)."""
        if center_frac is None:
            center_frac = self.sleeve_center_frac
        r_mid = self.r_pipe + center_frac * (self.r_wall - self.r_pipe)
        return r_mid - self.t_sleeve / 2, r_mid + self.t_sleeve / 2
