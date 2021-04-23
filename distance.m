function [out cost_flag] = distance(in)
% Distance criterion for estimation
global pdfdata adjtype          %#ok<NUSED>

lbar       = in(1);
alpha      = in(2);
ksi        = in(3);
rho        = in(4);
stdMC      = in(5);

version      = 3;
gridSpread   = 0.1;     % extra spread for price grid as a share of (PMAX-PMIN)
finegrid     = 0;       % fine (1) or coarse (0) grid
accuracy     = 0;       % accuracy of SS calculation from 0 (lowest) to 4 (highest). Controls SS error tolerances.

p_iter;


% this gives to freqpchanges the same weight as that of the entire vector
% of histogram counts (or 25 times more weight than any single count)
 out = length(prob)*norm(freqpchanges-0.10) + norm(prob-pdfdata);        

cost_flag = 1;



% out = length(prob)*norm(freqpchanges-0.205) + norm(prob-pdfdata);           %FREQ 20.5pct  (MidAC)
% out = 1*norm(freqpchanges-0.205) + norm(prob-pdfdata);           %FREQ 20.5pct  (MidAC)

% alternatively, one may want to minimize the max distance
% out = max(abs(freqpchanges-0.10) ; abs(prob-nielN));   % infinity norm

% this gives to freqpchanges the same weight as any single histogram count
% out = (freqpchanges-0.10)^2 + (prob-nielN).^2      ;   % least squares