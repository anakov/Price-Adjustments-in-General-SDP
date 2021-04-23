% Computes and plots impulse-response functions
clear STATEHISTORY; clear JumpHistory 
disp(sprintf('\n'))  
TT = 20;              %periods 1:19 fit nicely in graph 0:20
INITCONDIT=0;

% SET MONEY IMPULSES
time1moneyShock = 0;
time2moneyShock = 10*jacstep; % shock should be approx size of jacstep to preserve linearity in aggregate shocks
time3moneyShock = 0;

% SPECIFY MONEY SHOCK PROCESS for periods 1:TT
moneyshocks = [time1moneyShock time2moneyShock time3moneyShock zeros(1,TT-3)];  % 3 shocks + zeros
shocktime = 2;
scalefactor = abs(1/moneyshocks(shocktime ));  % for plots expressed in pp deviations (not % difference from SS)

distsim_jmcb;
compute_IRFs_jmcb;
plot_IRFs_jmcb;

% checkplot;
% showmovie;
