function R = ccmSlipModel(P, b)
%CCMSLIPMODEL  Model A: 0-D lumped force-balance / lubrication model
%for slip-enhanced close-contact melting (sCCM) of a solid PCM column
%sinking under gravity around a vertical, pulsed cylindrical heater
%sleeve.
%
%   R = CCMSLIPMODEL(P, b) runs the model for slip length b [m] with
%   parameters P (see parameters.m) and returns a struct R with fields
%   t, Hs, meltFraction, delta, qH, vSink, pulseOn, capped, milestones,
%   b.
%
%PHYSICAL PICTURE
%-----------------
%The tube is vertical, so unlike the horizontal grooved-plate CCM setup
%(gravity normal to the melting plate), here gravity acts *along* the
%tube axis. The heater is a sleeve running up the axis; solid PCM fills
%the annulus around it and, over the whole column height, well above
%the heated section too. As PCM melts off the sleeve's surface, a thin
%lubricating melt film separates the sleeve from the (still solid,
%still rigid) PCM column above it. The column's own net weight (gravity
%minus buoyancy) presses down on this film and squeezes melt axially
%downward out of the heated section; mass conservation (melting flux in
%= axial drainage flux out) plus a lubrication force balance sets the
%film's quasi-steady thickness delta(t), from which the heat flux
%q = kLiq*dT/delta follows. This is the vertical-tube analogue of
%classic horizontal CCM (Bejan; Moallemi & Viskanta 1985) and reduces
%to it algebraically -- see "Sanity check" below.
%
%The slip coating (Li et al. 2026, Nature, "Pulse heating and slip
%enhance charging of phase-change thermal batteries") does not change
%the conduction resistance across the film directly (heat still crosses
%by conduction, q = kLiq*dT/delta, exactly as without slip) -- it
%changes the *equilibrium* delta the force balance settles on, by
%reducing the viscous resistance to the axial drainage flow that must
%be squeezed out through the film. A thinner equilibrium film means a
%higher q. This matches the mechanism reported in the reference paper.
%
%DERIVATION (plane-unrolled lubrication, one slip wall)
%-------------------------------------------------------
%Unroll the annular film (width delta << rH) into a 2-D channel of gap
%delta, spanwise width w = 2*pi*rH, running the wetted height Hc. Plane
%Poiseuille flow driven by -dp/dz = G, with Navier slip length b on the
%heater wall (y=0) and no-slip on the solid's melting surface (y=delta):
%
%   u(y) = -(G/2mu) y^2 + C1 y + C2,   C1 = G delta^2 / (2 mu (delta+b)),
%   C2 = b C1
%
%Integrating u over the gap gives the flow rate per unit width
%
%   Q' = G * delta^3 (delta + 4b) / (12 mu (delta + b))
%      = G * delta_eff^3 / (12 mu),   delta_eff^3 := delta^3(delta+4b)/(delta+b)
%
%which recovers the ordinary no-slip result Q' = G delta^3/(12 mu) at
%b=0, and the b -> infinity (shear-free wall) limit Q' = G delta^3/(3 mu)
%(4x enhancement), both standard checks on a one-wall-slip plane
%Poiseuille solution (see deltaEffCubed.m).
%
%Melting is distributed uniformly over the wetted height Hc (uniform-dT,
%uniform-delta lumped approximation), so the axial volumetric flow at
%height z (measured from the open bottom, z=0, up to the closed top of
%the wetted zone, z=Hc) is the melt generated above it:
%
%   Q(z) = q_h w (Hc - z) / (rhoLiq*Lf),   q_h = kLiq*dT/delta
%
%Using dp/dz = 12*mu*Q(z) / (w*delta_eff^3) and integrating twice gives
%the excess pressure profile p(z), and integrating p(z) over the
%contact area gives the total upward force supporting the column's net
%weight Wnet:
%
%   Wnet = 4*muLiq*w*q_h*Hc^3 / (delta_eff^3 * rhoLiq*Lf)             (*)
%
%With q_h = kLiq*dT/delta substituted, (*) is one nonlinear equation in
%the one unknown delta (given Wnet, dT, Hc, b); solved numerically per
%time step by solveDelta.m.
%
%Sanity check (b=0):  delta * delta_eff^3 = delta^4, so
%   delta = [ 4*muLiq*w*kLiq*dT*Hc^3 / (Wnet*rhoLiq*Lf) ]^(1/4)
%which is exactly the classic CCM 1/4-power scaling for film thickness
%(propto (mu*k*dT*L^3 / (W*rho*Lf))^(1/4)) reported throughout the CCM
%literature (Bejan 1994; Moallemi & Viskanta 1985) -- confirms the
%lumped derivation above is dimensionally and structurally consistent
%with the established (no-slip) theory before slip is added.
%
%MODEL ASSUMPTIONS / LIMITATIONS (stated explicitly, not hidden)
%------------------------------------------------------------------
%  * Lumped in z: dT and delta are treated as spatially uniform over
%    the wetted height Hc at each instant (only the *flow rate* Q(z) is
%    allowed to vary with z, which is what the force balance needs). A
%    fully-resolved sCCM model would let delta vary with z; this
%    lumped version is a standard first simplification in the CCM
%    literature and is adequate for overall charge-time trends.
%  * Pulsing: the thin film's own thermal mass is neglected, so during
%    an OFF interval (wall adiabatic) dT -> 0 instantly and melting
%    pauses (q_h=0) rather than decaying gradually. This is
%    conservative in the sense that it does not credit the model with
%    any "coasting" melting during OFF.
%  * The force balance is algebraic (quasi-steady), not an ODE for
%    delta itself, so delta(t=0) is simply whatever solves (*) at the
%    initial (finite) Wnet -- there is no zero-gap start-up singularity
%    to regularize. What the model *does* skip over is the brief real
%    conduction-limited transient before quasi-steady lubrication is
%    established (well documented in the CCM literature, e.g. Moallemi
%    & Viskanta 1985); this is short compared to the charge times of
%    interest here and is not separately modeled. P.delta0 instead
%    serves as a minimum achievable film thickness (surface roughness /
%    non-ideal-contact floor): delta is clamped to it, which in turn
%    caps q_h wherever the smooth-wall force balance would otherwise
%    want an unrealistically thin film.
%  * Breaks down (delta -> unbounded) as Wnet -> 0 near complete
%    melt-out (Hs -> 0); capped at solveDelta's upper bound and flagged
%    in R.capped. t99 is reported (see milestonesFromSeries.m) but
%    flagged unreliable for this model -- t95 is the more robust
%    summary number.

w  = 2*pi*P.rH;
Ac = P.Across;
dt = P.dtA;
nSteps  = floor(P.tMax/dt) + 1;
HsFloor = 1e-4;   % treat as "fully melted" below this remaining height

t       = zeros(nSteps,1);
Hs      = zeros(nSteps,1);
delta   = nan(nSteps,1);
qH      = zeros(nSteps,1);
vSink   = zeros(nSteps,1);
pulseOn = false(nSteps,1);
capped  = false(nSteps,1);

HsCur = P.Hcol;
nUsed = nSteps;

for i = 1:nSteps
    ti = (i-1)*dt;   % computed fresh each step (not accumulated) to
                     % avoid float drift landing on the wrong side of a
                     % pulse edge
    t(i)  = ti;
    Hs(i) = HsCur;

    isOn = pulseState(ti, P);
    pulseOn(i) = isOn;
    Hc = min(HsCur, P.Hheater);

    if isOn && Hc > 0 && HsCur > HsFloor
        dT   = P.Th - P.Tm;
        Wnet = (P.rhoSol - P.rhoLiq) * P.g * Ac * HsCur;
        if Wnet > 0
            d = solveDelta(Wnet, P.muLiq, P.kLiq, dT, Hc, w, P.rhoLiq, P.Lf, b);
            capped(i) = (d <= 1e-9*1.0001) || (d >= 5e-2*0.9999);
            d = max(d, P.delta0);   % surface-roughness / non-ideal-contact floor
            q = P.kLiq * dT / d;
        else
            d = NaN; q = 0;
        end
    else
        d = NaN; q = 0;
    end

    delta(i) = d;
    qH(i)    = q;

    if Hc > 0
        v = q * (2*pi*P.rH*Hc) / (P.rhoSol*P.Lf*Ac);
    else
        v = 0;
    end
    vSink(i) = v;

    if HsCur <= HsFloor
        nUsed = i;
        break;
    end

    HsCur = max(HsCur - v*dt, 0);
end

t       = t(1:nUsed);
Hs      = Hs(1:nUsed);
delta   = delta(1:nUsed);
qH      = qH(1:nUsed);
vSink   = vSink(1:nUsed);
pulseOn = pulseOn(1:nUsed);
capped  = capped(1:nUsed);

meltFraction = 1 - Hs/P.Hcol;
milestones   = milestonesFromSeries(t, meltFraction);

R.t            = t;
R.Hs           = Hs;
R.meltFraction = meltFraction;
R.delta        = delta;
R.qH           = qH;
R.vSink        = vSink;
R.pulseOn      = pulseOn;
R.capped       = capped;
R.milestones   = milestones;
R.b            = b;

end
