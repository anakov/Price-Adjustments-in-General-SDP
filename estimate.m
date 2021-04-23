% Estimates parameters of adjustment function and idiosyncratic shocks
clear; clc
global pdfdata adjtype         
set(0,'DefaultFigureWindowStyle','normal')  % undock all figures
load acnielsen

lbound = -0.5;
hbound =  0.5; 
edges  = [-inf linspace(lbound,hbound,24) inf];   % same bin edges as in calcstats; produces 25 bins

pdfdata = histc(data,edges);
pdfdata = pdfdata(1:end-1);    % the last bin counts any values that match EDGES(end); see help histc
pdfdata = pdfdata./sum(pdfdata);

% lambda(L) = lbar / (1-lbar)*(alpha/L)^ksi + lbar       where L = D/wbar;
% lbar  = NaN   ;  
% alpha = 0.03;  % alpha: when L=alpha, hazard(d)=lambdabar (like "menu cost")
% ksi   = NaN   ;  % ksi=0 --> Calvo, ksi=inf --> Menu cost   (like "info cost" in Woodford 2008)
% rho   = 0.87;  % rho
% stdMC = 0.12;  % standard deviation of marginal cost
   
lbar  = 0.108928579454875;
alpha = 0.031083222345823;
ksi   = 0.293675662699625;
rho   = 0.881156092163841;
stdMC = 0.147421080785644;

lb    = [0.10 0.02 0.2 0.85  0.13]';
ub    = [0.12 0.04 0.4 0.93  0.18]';

adjtype = 0;        % 0-SSDP; 2-Calvo; 3-GL; 4-Woodford; 5-DKW 

guess = [lbar alpha, ksi, rho, stdMC]';   
options=optimset('Display','iter','OutputFcn', @optimprint,...
                  'TolFun',1e-5,'TolX',1e-5,'DiffMaxChange',5e-1,'DiffMinChange',1e-7);
 [Paramvector,fval,exitflag,output,lambda,grad] = fmincon('distance',guess,[],[],[],[],lb,ub,[],options);

% Alternative procedure: simulated annealing
% asamin('set','rand_seed',696969);
% asamin('set','asa_out_file','estimate.log');
% asamin('set','test_in_cost_func',0);
% [fstar,xstar,grad,hessian,state] = asamin ('minimize', 'distance', guess, lb, ub, -ones(size(guess)))

