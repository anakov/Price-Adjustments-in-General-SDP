% Solves value function by iteration on a grid with interpolation
% version 17 april 2008, written for model detrended by M

VDIFF=inf;                                                 % reset difference in value function
Viter=0;                                                   % reset counter

while (VDIFF>Vtol*max([1 10^(3-piter)]) && Viter<10000)    % iterate to convergence of V 
    Viter=Viter+1;                                         % increment counter
    if rem(Viter,50)==0  && showconverge                   % report convergence progress
        convergence_report(adjtype,finegrid,MMiter,...
            sizeofMMgrid,pDIFF,piter,VDIFF,Viter,[],[]); 
    end
    
    % Time t+1 values
              
    %  Quadratic value function interpolation at each state
    if RECURSE==0
       % Discrete approximation for first N iterations (to speed up the convergence)
       if (Viter<=800)                                       
          [M,OPTindex] = max(V);                               % MAX of discrete approximation
          pstar = Pgrid(OPTindex)';                            % M = max V; OPTindex is position of best policy
       end
       for col=1:nums                                          % Loop over states (columns of V)
          [maxV, maxind] = max(V(:,col));                      % Discrete maximum and index of maximant
          if nums==1 && mu==1                                  % if rep agent, zero SS inflation
              M = maxV;                                        % no interpolation: opt price is flex price
              pstar = Pgrid(maxind);
          elseif maxind==1 || maxind==nump
              M(col) = maxV;
              pstar = Pgrid(maxind);
          else
              localx = Pgrid(maxind-1:maxind+1);               % Local grid of 3 points surrounding maximant
              localV = V(maxind-1:maxind+1,col);               % Value of V on local grid
              XMAT = [ones(3,1) localx localx.^2];             % Build regressors
              betaquad = (XMAT'*XMAT)\XMAT'*localV;            % OLS quadratic fit to V over local grid
              pstar(col) = -0.5*betaquad(2)/betaquad(3);       % Maximant of interpolant
              M(col) = [1 pstar(col) pstar(col)^2]*betaquad;   % Maximum of interpolant
          end
       end
    elseif RECURSE==1
    %  error('ERROR! interpolation program not valid for case of (p,pi) plans!')
       [M,OPTindex] = max(V);                              % MAX of discrete approximation
       pstar = Plongvec(OPTindex)';                        % M = max V; OPTindex is position of best policy
       pistar = pilongvec(OPTindex);
    end

    
    D = ones(numlong,1)*M - V;                             % D is the value of adjustment
    if any(any(D<-eps^.5)),
       % disp('ERROR: Negative D after quadratic interpolation of V. Stretch out grid')
    end
    D(D<0) = 0;
    
    [lambda, POSSIBILITIES] = adjustment(adjtype, V, D, wbar, ksi, alpha, lbar);

  % Time t values
    iterV = PAYOFFMAT + beta*RECURSEMAT'*POSSIBILITIES*TRANSMAT;   
                                                           % iterV is current payoff plus disc. continuation value
    VDIFF = max(max(abs(iterV-V)));                        % change in value function (sup norm)
    V = iterV;                                             % updating V

end

if ( any(pstar == PMIN) || any(pstar == PMAX) )
    disp('Maximum attained at grid boundary. May need to stretch out grid')
end


