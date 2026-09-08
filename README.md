# simtry
CCM Simulation

| Folder | What it is |
|---|---|
| [`PCM_MATLAB/`](PCM_MATLAB/) | Version 1: 3D conduction + phase-change MATLAB/Octave solver for a PCM store heated by two concentric copper heater sleeves (no gravity-driven melt motion). |
| [`PCM_sCCM_MATLAB/`](PCM_sCCM_MATLAB/) | Version 2 (MATLAB/Octave): same sCCM models as below, ported to MATLAB in `PCM_MATLAB`'s per-function-per-file style. **Use this one if you work in MATLAB.** |
| [`PCM_sCCM/`](PCM_sCCM/) | Version 2 (Python): models of slip-enhanced close-contact melting (sCCM) in a vertical pulsed-heater tube — a 0-D sinking-column force-balance model and a 1-D radial transient enthalpy PDE, grounded on Li et al. 2026 (*Nature*). |
