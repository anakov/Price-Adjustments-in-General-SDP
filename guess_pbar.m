% Initial guess for the price level 
% p=P/M
switch   mu
    case muGLhigh
             pbar = 0.2882;
             Calvoparam = 0.205; 
    case muGanHi
             pbar = 0.2882;
             Calvoparam = 0.25; 
    case muDomin
             pbar = 0.0223;                                     
             Calvoparam = 0.10;
    case muACNiel
             pbar = 0.0223;                                     
%             Calvoparam = 0.10;
             Calvoparam = 0.205;
    otherwise
             pbar = wbar*epsilon/(epsilon-1);  % average flex price markup over average marginal cost
             Calvoparam = 0.10; 
end

% pbar = 5*wbar; % check for uniqueness with different intial values