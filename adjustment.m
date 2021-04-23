function [lambda, POSSIBILITIES] = adjustment(adjtype, V, D, wage, ksi, alpha, lbar)
% Computes adjustment probability as a function of the gain from adjustment

L =  D./wage;                                          % Adjustment depends on D/w, value in units of labor time
switch adjtype
    case 0                                             % BASELINE SSDP MODEL 
      aL = (alpha./L).^ksi;                                % 1 when ksi=0; 0 or inf when ksi=inf
      lambda = lbar./((1-lbar)*aL+lbar);                   % probability of adjustment as a function of the state 
    case 1                                             % SIMPLER SSDP MODEL
      Lksi = L.^ksi;                                       % L raised to ksi 
      lambda = Lksi./(lbar+Lksi);                          % probability of adjustment as a function of the state 
    case 2                                             % CALVO MODEL 
      lambda = lbar*ones(size(D));                     %   Constant probability of adjustment
    case 3                                             % GOLOSOV-LUCAS MODEL
      lambda = lamcontin(L,alpha);                         % probability of adjustment: 0 or 1
    case 4                                             % WOODFORD'S MODEL
      argexp=(alpha-L)*ksi;                                % argument in exponent function
      lambda = lbar./((1-lbar).*exp(argexp) + lbar);       % probability of adjustment as a function of the state 
    case 5                                             % DOTSEY-KING-WOLMAN MODEL
      aL = (alpha./L).^ksi;                                % 1 when ksi=0; 0 or inf when ksi=inf
      lambda = lbar./((1-lbar)*aL+lbar);                   % probability of adjustment as a function of the state 
    case 6                                             % TRUNCATED CALVO 
      lambda = lamcontin(L,alpha);                         % probability of adjustment 1 outside bands
      lambda(lambda==0) = lbar;                            % constant positive probability inside bands
    case 7
      lambda = (lbar+L)./(1+L);                        % SIMPLE SSDP with lambda(0)>0
    case 8                                             % BASELINE SSDP with lambda(0)>0 
      lambda = (lam0 + lbar*L.^ksi)./...
          ((1-lbar)*alpha^ksi + lbar*L.^ksi);  
end

if nargout==2
   if adjtype<3 || adjtype>5
     alpha = 0;                                        % No fixed cost is subtracted in Calvo and SSDP models
   elseif adjtype==5
     alpha = expectmc(adjtype, D, wage, ksi, alpha, lbar); % Expected menu cost in DKW model 
   end
   POSSIBILITIES = V + (D - alpha*wage).*lambda;       % Continuation value in all models
end
