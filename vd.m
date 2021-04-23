clear STATEHISTORY; clear JumpHistory 
% Computes variance decomposition and Phillips curve regression

fprintf('\n')  
TT = 3*100;      % sample size (months; must divide by 3)
INITCONDIT=0;

% SPECIFY MONEY SHOCK PROCESS for periods 1:TT
randn('state',0)
shocksize = jacstep;
scalefactor = 1/shocksize;
moneyshocks = [shocksize *randn(1,TT)];  %#ok<NBRAK> %simulating random history
time1moneyShock = moneyshocks (1);

TFPshocks = zeros(1,TT);
time1TFPshock = TFPshocks(1);

distsim_jmcb;
compute_IRFs_jmcb;    % run money growth shock

if rem(TT,3)==0 
    
% convert to quarterly frequency
  C_pathQ = mean(reshape(C_path,3,TT/3));
  PI_identity_meanQ = mean(reshape(PI_identity_mean+mu,3,TT/3));
  d_R_pathQ = mean(reshape(d_mu_path,3,TT/3));
  
  pchCQ = C_pathQ(2:end)./C_pathQ(1:end-1)-1;

% gdp growth and delfator inflation quarterly s.d. during 1984-2008
  data_V_cons_growth = 0.00510;  % okay:  cgg_mpr_qje.wf1
  data_V_infl        = 0.00246;  % okay:  cgg_mpr_qje.wf1 

  scaleUp = data_V_infl/std(PI_identity_meanQ);          % scale up to explain all observed inflation
  
  std_PI_identity_meanQ  = std(PI_identity_meanQ)*scaleUp;
  std_pchCQ     = std(pchCQ)*scaleUp;

  VD_inflation   = std_PI_identity_meanQ/data_V_infl  ;  % equals 1 by construction
  VD_cons_growth = std_pchCQ/data_V_cons_growth ;  
  
% 2SLS regression of consumption on inflation     
% first stage regression: inflation on exogenous shock 
  regressors = [ones(size(d_R_pathQ')) d_R_pathQ'];
  if rank(regressors)==size(regressors,2);
     B = regress(PI_identity_meanQ',regressors);
  else
     error('Colinear regressors');
  end
  PI_projQ = B(1) + B(2)*regressors(:,2);
  
% second stage regression: output on predicted inflation     
  regressors = [ones(size(PI_projQ)) 4*log(PI_projQ)];
  if rank(regressors)==size(regressors,2);
     [B,BINT,R,RINT,STATS] = regress(log(C_pathQ'),regressors );
  else
     error('Colinear regressors');
  end

  
% Print output
  fprintf('100 x std dev of money shock: %0.3g',100*std(moneyshocks)*scaleUp)
  fprintf('\n')
  fprintf('Model implied 100 x std of quarterly inflation             : %0.3g \n',100*std_PI_identity_meanQ)
  fprintf('Actual 100 x std of quarterly deflator inflation 1984-2008 : %0.3g \n',100*data_V_infl)
  fprintf('Share of inflation variance due to money growth shocks     : %0.3g%% \n', 100*VD_inflation)
  fprintf('\n')
  fprintf('Model implied 100 x std of quarterly output growth         : %0.3g \n',100*std_pchCQ)
  fprintf('Actual 100 x std of quarterly output growth 1984-2008      : %0.3g \n',100*data_V_cons_growth)
  fprintf('Share of output variance due to money growth shocks        : %0.3g%% \n', 100*VD_cons_growth)
  fprintf('\n')
  fprintf('Phillips curve: log(C_pathQ) = alpha + beta(4log(PI_projQ)) \n')
  fprintf('Data frequency: quarterly (average of monthly simulated observations) \n')
  fprintf('Estimation method: 2SLS \n')
  fprintf('Instrument for inflation: exogenous aggregate shock \n')
  fprintf('\n')
  fprintf('Slope coefficient beta                                     : %0.3g \n', B(2))
  fprintf('Standard error for slope coefficient                       : %0.3g \n', abs(B(2)-BINT(2,1))/2)
  fprintf('R2 of regression                                           : %0.3g \n', STATS(1))
  fprintf('\n')
end
