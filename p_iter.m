% Computes general equilibrium steady-state on a grid
% version 17 april 2008, model detrended by M; iteration on pbar
tic
param;                                                           % set parameters including initial guess for pbar
GE_setmat;                                                       % construct initial matrices
pDIFF=inf;                                                       % reset difference in p
piter=0;                                                         % reset counter
  
while (pDIFF>ptol && piter<20)                                   % iterate to convergence of p
    piter=piter+1;                                               % increment counter
    Cbar = ((1-beta/mu)/(nu*pbar))^(1/gamma);                    % calculate Cbar from pbar
    PAYOFFMAT = Cbar*(pbar^epsilon)*...
               (PMAT.^(1-epsilon)-wbar*sMAT.*PMAT.^(-epsilon));  % compute payoff matrix
    GE_V_fun_iter;                                               % get value function given p, C, w
    GE_Pdist_iter;                                               % get stationary distribution (p,mc) given C
    newp = sum(exp(Plongvec).^(1-epsilon).*sum(Pdist,2)).^(1/(1-epsilon)) ;   % compute new p     
    pDIFF=abs(newp-pbar);                                        % compute norm
    pbar=newp;                                                   % update p
end

% PRINT OUTPUT 
calcstats                                                    % compute statistics
%  printstats                                                   % print out statistics
%  plotfigs                                                     % plot figures

if pDIFF>ptol, disp(sprintf('WARNING: Convergence failed for p. pDIFF=%d. Stretch out grid?',pDIFF)), end
if VDIFF>Vtol, disp(sprintf('WARNING: Convergence failed for V. VDIFF=%d. Stretch out grid?',VDIFF)), end
if PdistDIFF>PdistDIFFtol, 
   disp(sprintf('WARNING: Convergence failed for Pdist. PdistDIFF=%d. Stretch out grid?',PdistDIFF)), end

