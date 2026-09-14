function rows = runPositionSweepRadial(P, matp, fracs, recordEvery, verbose)
%RUNPOSITIONSWEEPRADIAL  Vary the copper sleeve's radial position across
%the pipe->wall gap (holding the "lowest power to reach Tm within the
%ramp target" control rule fixed at each position) to find the
%placement that minimizes total melting time.
%
%   rows = RUNPOSITIONSWEEPRADIAL(P, matp, fracs, recordEvery, verbose)
%
%   Returns an array of structs, one per fraction tried, each with
%   frac, rIn, rOut, rMid, Pmin, t50, t90, t99, t999, dutyCycle.

if nargin < 3 || isempty(fracs), fracs = linspace(0.15, 0.85, 8); end
if nargin < 4 || isempty(recordEvery), recordEvery = 60; end
if nargin < 5 || isempty(verbose), verbose = true; end

rows = struct('frac', {}, 'rIn', {}, 'rOut', {}, 'rMid', {}, 'Pmin', {}, ...
              't50', {}, 't90', {}, 't99', {}, 't999', {}, 'dutyCycle', {}, 'converged', {});

for i = 1:numel(fracs)
    frac = fracs(i);
    G = createGeometryRadial(P, frac);
    [Pmin, ~] = findMinPowerRadial(P, G, matp);
    R = runFullMeltRadial(P, G, matp, Pmin, P.Tm, recordEvery, false);

    row.frac = frac;
    row.rIn = G.rSleeveIn; row.rOut = G.rSleeveOut;
    row.rMid = 0.5*(G.rSleeveIn + G.rSleeveOut);
    row.Pmin = Pmin;
    row.t50 = R.t50; row.t90 = R.t90; row.t99 = R.t99; row.t999 = R.t999;
    row.dutyCycle = R.dutyCycle;
    row.converged = R.converged;
    rows(end+1) = row; %#ok<AGROW>

    if verbose
        printf(['  frac=%4.2f  rMid=%6.2fmm  Pmin=%7.2f W/m  ', ...
                't50=%8s  t90=%8s  t99=%8s  t999=%8s\n'], ...
            frac, row.rMid*1000, Pmin, ...
            numOrNA(R.t50), numOrNA(R.t90), numOrNA(R.t99), numOrNA(R.t999));
    end
end

end

function s = numOrNA(x)
if isnan(x)
    s = 'n/a';
else
    s = sprintf('%.1f', x);
end
end
