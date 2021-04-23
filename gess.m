% Main steady-state program: computes general equilibrium steady-state and reports statistics
% Programs for "Price Adjustments in a General Model of
% State-Dependent Pricing", Working paper 0824, Bank of Spain
% Copyright (C) James Costain and Anton Nakov (2008)

% PROGRAM EXECUTION PARAMETERS                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
version      = 1;                              % Program version (1-all estimated; 2-fixed productivity; 3-estimation run)
                                               %    Set 3 if estimating the model with estimate.m. Otherwise set 1 to
                                               %    simulate an estimate in which productivity has been left free,
                                               %    or 2 to simulate an estimate in which productivity was fixed.
if version<3, adjtype = 0; end                 % 0-SSDP; 2-Calvo; 3-Fixed menu cost; 4-Woodford; 5-Stochastic menu costs
gridSpread   = 0.1;                            % extra spread for price grid as a share of (PMAX-PMIN)
finegrid     = 1*(adjtype==0);                 % fine (1) or coarse (0) grid
accuracy     = 0;                              % accuracy of SS calculation from 0 (lowest) to 4 (highest). Controls SS error tolerances.

p_iter
printstats
plotfigs