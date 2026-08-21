"""
Driver script: reproduce paper-style CCM-on-grooved-surfaces figures.

Runs the slip-modified CCM model (see model.py for the physics and its
documented approximations relative to the source paper) for:

  - a smooth (no-slip) reference surface,
  - three grooved surfaces spanning the paper's reported regimes
    (enhancing, negligible, suppressing melting relative to smooth),

at two Stefan numbers, each with the meniscus-curvature feedback ON
("refined model") and OFF ("conventional / constant-slip model"), and
produces:

  figures/melting_fraction.png  -- H/H0 vs normalized time   (like paper Fig. 3a-b)
  figures/heat_flux.png         -- average heat flux vs time (like paper Fig. 3c-d)
  figures/summary.txt           -- melting-rate enhancement/suppression table

Run:  python3 -m ccm_slip.run_paper_like
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ccm_slip import model as m

OUTDIR = os.path.join(os.path.dirname(__file__), "figures")
os.makedirs(OUTDIR, exist_ok=True)

L = 0.02       # m, block length along the drainage direction (paper: 20 mm)
H0 = 0.06      # m, initial solid height (paper: 60 mm)
H_STOP = 0.03  # stop integration once H falls to 3% of H0

SMOOTH = m.Groove(pitch=1.0, phi=0.0, name="Smooth (no-slip)")

# Showcase surfaces, chosen (via pitch_for_lambda0) to land in the three
# qualitative regimes this model actually predicts -- mild enhancement,
# near-neutral, and strong suppression -- relative to the no-slip film
# thickness h0 ~= 44-126 um spanned over the run (see README for the
# regime analysis). These are NOT the paper's own Surface A/B/C (whose
# exact slip-length closure is unavailable to this implementation); they
# demonstrate the same qualitative story -- enhance / negligible / suppress
# -- that the paper reports, using this model's own physics.
SURFACES = [
    m.Groove(pitch=m.pitch_for_lambda0(12e-6, 0.3), phi=0.3,
             name="Enhancing (l=0.33 mm, phi=0.3)"),
    m.Groove(pitch=m.pitch_for_lambda0(45e-6, 0.3), phi=0.3,
             name="Moderate (l=1.23 mm, phi=0.3)"),
    m.Groove(pitch=m.pitch_for_lambda0(140e-6, 0.5), phi=0.5,
             name="Suppressing (l=1.27 mm, phi=0.5)"),
]

ST_LIST = [0.026, 0.051]

COLORS = {"Smooth (no-slip)": "0.35",
          "Enhancing (l=0.33 mm, phi=0.3)": "#d62728",
          "Moderate (l=1.23 mm, phi=0.3)": "#2ca02c",
          "Suppressing (l=1.27 mm, phi=0.5)": "#1f77b4"}


def run_all():
    results = {}
    for St in ST_LIST:
        dT = m.dT_from_stefan(St)
        for surf in [SMOOTH] + SURFACES:
            for meniscus in (False, True):
                key = (St, surf.name, meniscus)
                print(f"running St={St:.3f}  {surf.name:32s}  "
                      f"meniscus={'on ' if meniscus else 'off'} ...", flush=True)
                res = m.run_ccm(surf, dT, L=L, H0=H0, meniscus=meniscus,
                                 H_stop_frac=H_STOP)
                results[key] = res
    return results


def melt_time_to_frac(res: m.SimResult, frac=0.5):
    """Interpolate the time at which H/H0 first drops to `frac`."""
    y = res.H / res.H0
    if y[-1] > frac:
        return np.nan
    idx = np.searchsorted(-y, -frac)
    if idx == 0:
        return res.t[0]
    t0, t1 = res.t[idx - 1], res.t[idx]
    y0, y1 = y[idx - 1], y[idx]
    return t0 + (y0 - frac) / (y0 - y1) * (t1 - t0)


def plot_melting_fraction(results):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.4), sharey=True)
    for ax, St in zip(axes, ST_LIST):
        dT = m.dT_from_stefan(St)
        t_ref = melt_time_to_frac(results[(St, SMOOTH.name, False)], frac=0.5)
        for surf in [SMOOTH] + SURFACES:
            for meniscus, ls in ((False, "--"), (True, "-")):
                res = results[(St, surf.name, meniscus)]
                tau = res.t / t_ref
                ax.plot(tau, res.H / res.H0, ls, color=COLORS[surf.name],
                         lw=1.8, alpha=0.9,
                         label=f"{surf.name}" if meniscus else None)
        ax.set_xlim(0, 2.0)
        ax.set_ylim(0, 1.02)
        ax.set_xlabel(r"$\tau = t / t_{ref}$  (smooth-surface half-melt time)")
        ax.set_title(f"St = {St:.3f}  ($\\Delta T$ = {dT:.2f} K)")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel(r"$H/H_0$")
    handles, labels = axes[0].get_legend_handles_labels()
    solid = plt.Line2D([], [], color="k", ls="-", label="refined (meniscus)")
    dashed = plt.Line2D([], [], color="k", ls="--", label="conventional (const. slip)")
    fig.legend(handles + [solid, dashed], labels + ["refined (meniscus)", "conventional (const. slip)"],
                loc="lower center", ncol=3, fontsize=8, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("CCM melting fraction on grooved slip surfaces\n"
                 "(own reduced-model reproduction -- see model.py docstring for scope/approximations)",
                 fontsize=10)
    fig.tight_layout(rect=[0, 0.13, 1, 0.92])
    path = os.path.join(OUTDIR, "melting_fraction.png")
    fig.savefig(path, dpi=160)
    print("wrote", path)


def plot_heat_flux(results):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    for ax, St in zip(axes, ST_LIST):
        t_ref = melt_time_to_frac(results[(St, SMOOTH.name, False)], frac=0.5)
        q_ref = results[(St, SMOOTH.name, False)].q_avg[0]
        for surf in [SMOOTH] + SURFACES:
            for meniscus, ls in ((False, "--"), (True, "-")):
                res = results[(St, surf.name, meniscus)]
                tau = res.t / t_ref
                ax.plot(tau, res.q_avg / q_ref, ls, color=COLORS[surf.name],
                         lw=1.8, alpha=0.9)
        ax.set_xlim(0, 1.0)
        ax.set_xlabel(r"$\tau = t / t_{ref}$")
        ax.set_title(f"St = {St:.3f}")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel(r"$q''/q''_{ref}$  (normalized by smooth-surface initial flux)")
    fig.suptitle("Average film heat flux on grooved slip surfaces\n"
                 "(own reduced-model reproduction -- see model.py docstring for scope/approximations)",
                 fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    path = os.path.join(OUTDIR, "heat_flux.png")
    fig.savefig(path, dpi=160)
    print("wrote", path)


def write_summary(results):
    lines = []
    lines.append("Melting-rate enhancement/suppression relative to the smooth surface")
    lines.append("(conventional model = constant slip length; refined model = with")
    lines.append(" meniscus-curvature feedback -- see model.py for what this means)")
    lines.append("")
    lines.append("For reference, the paper (whose exact slip-length closure is not")
    lines.append("available to this model -- see model.py) reports, across its own")
    lines.append("three surfaces at the same superheat: melting enhanced by >15% on")
    lines.append("one grooved surface, suppressed by ~21% on another, and negligibly")
    lines.append("changed on a third. The three surfaces below are chosen (by this")
    lines.append("model's own slip-vs-film-thickness ratio) to demonstrate the same")
    lines.append("three qualitative regimes, not to reproduce those exact numbers.")
    lines.append("")
    for St in ST_LIST:
        lines.append(f"St = {St:.3f}")
        t_smooth = melt_time_to_frac(results[(St, SMOOTH.name, False)], frac=0.5)
        for surf in SURFACES:
            for meniscus in (False, True):
                res = results[(St, surf.name, meniscus)]
                t50 = melt_time_to_frac(res, frac=0.5)
                pct = (t_smooth / t50 - 1.0) * 100.0  # faster half-melt = enhanced rate
                tag = "refined " if meniscus else "convent."
                lines.append(f"  [{tag}] {surf.name:32s} half-melt time change: {pct:+6.1f}%")
        lines.append("")
    text = "\n".join(lines)
    path = os.path.join(OUTDIR, "summary.txt")
    with open(path, "w") as f:
        f.write(text + "\n")
    print(text)
    print("wrote", path)


if __name__ == "__main__":
    results = run_all()
    plot_melting_fraction(results)
    plot_heat_flux(results)
    write_summary(results)
