
MATLAB Program files for the manuscript "Price Adjustments in a General Model of State-Dependent Pricing"
Bank of Spain working paper 0824, November 2008 
(C) James Costain and Anton Nakov


MAIN PROGRAMS TO RUN

  gess               - Main steady-state program: computes general equilibrium steady-state and reports statistics
  gedyn              - Main dynamics program: computes general equilibrium dynamics and shows impulse-responses
  estimate           - Estimates parameters of adjustment function and idiosyncratic shocks


COMPUTATIONAL DETAILS

Program gess.m (steady-state calculation) takes less than 1 minute on an ordinary Pentium 4 with 3Ghz CPU and 1GB of memory.
Program gedyn.m (dynamics solution) takes around 3-4 minutes on the same machine. 
The slowest step is the QZ decomposition, which uses MATLAB's internal function qz.m.


LIST OF ALL PROGRAMS IN ALPHABETICAL ORDER

  adjustment         - Computes adjustment probability as a function of the gain from adjustment
  calcstats          - Computes steady-state statistics 
  compute_IRFs_jmcb  - Computes dynamic paths as a function of initial state and shocks
  convergence_report - Reports convergence progress
  disclyap           - Solves discrete Lyapunov equation
  distance           - Distance criterion for estimation
  distsim_jmcb       - Computes impulse-responses based on Klein's state-space solution
  dyneqklein_jmcb    - Defines equations of stochastic model
  dynsolveklein_jmcb - Solves dynamic general equilibrium
  estimate           - Estimates parameters of adjustment function and idiosyncratic shocks
  expectmc           - Computes the expected menu cost in the stochastic menu cost model
  GE_Pdist_iter      - Computes the steady state distribution of firms on (price, productivity)
  GE_setmat          - Sets up matrices for general equilibrium model
  GE_V_fun_iter      - Solves value function by iteration on a grid with interpolation
  gedyn              - Main dynamics program: computes general equilibrium dynamics and shows impulse-responses
  gess               - Main steady-state program: computes general equilibrium steady-state and reports statistics
  guess_pbar         - Initial guess for the price level 
  hpfilter           - Hodrick-Prescott filter 
  infdecomp          - Decomposes inflation response into intensive, extensive, and selection effects
  irf_jmcb           - Computes and plots impulse-response functions
  jacob_reiter       - Computes Jacobian by forward differences
  kleinsolve_jmcb    - Implements Klein's QZ decomposition method for solving linear RE models
  ks                 - Computes Kolmogorov-Smirnov statistic for equality of two cdf's
  lamcontin          - Interpolates step adjustment function for the menu cost model
  midriplot          - Plots the histogram of price changes for the AC Nielsen dataset of Midrigan
  optimprint         - Prints current parameter vector at each iteration of the estimation procedure
  p_iter             - Computes general equilibrium steady-state on a grid
  param              - Sets model parameters and sets up grid
  plot_IRFs_jmcb     - Plots impulse-response functions
  plotfigs           - Plots steady-state distribution, histogram, value function and related objects
  printstats         - Prints out steady-state statistics
  tauchen            - Converts a VAR(1) into a Markov-Chain using Tauchen's method
  vd                 - Computes variance decomposition and Phillips curve regression
