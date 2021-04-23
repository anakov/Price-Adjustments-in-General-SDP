% Defines equations of stochastic model
% Input: X: All variables, current and lagged
% Outputs: Equation residuals
% WARNING: this setup is written assuming STICKY PRICES (ie RECURSE=0),
% not sticky policies

function Resid = dyneqklein(Y)
global Params;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD parameters
 
PMAT = Params.PMAT; 
sMAT = Params.sMAT; 
gamma = Params.gamma; 
chi = Params.chi; 
epsilon= Params.epsilon; 
adjtype = Params.adjtype; 
alpha = Params.alpha; 
lbar = Params.lbar; 
ksi = Params.ksi; 
beta = Params.beta; 
mu = Params.mu; 
nu = Params.nu; 
TRANSMAT = Params.TRANSMAT; 
Pgrid = Params.Pgrid; 
pstep = Params.pstep; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KLEIN setup: X=[V;c;p;Phihat] 
% z is the percent deviation of mu: 
%  mu_t = mu*exp(z_t) =approx= mu*(1 + z_t), where z_t+1 = phiz z_t + iidshock
% PhiHatNow is dist at beginning of t before firms adjust
%            AND before money shock z_t is realized.
% I.e. PhiHatNow_t = Phi_t-1 = end of period t-1 production distribution
%      
% State variables are: PhiHatNow and today's money process z_t

% LOAD VARIABLES 
  nV = Params.nV;
  nPhi = Params.nPhi;
% nz = Params.nz; % not used at present
  
% BREAK the Y vector into four parts:
  xnow = Y(Params.ix);         % variables now
  xnext = Y(Params.ixnext);    % variables next

  znow = Y(Params.iz);       % exogenous shock process (possibly correlated)
  znext = Y(Params.iznext);  % ShockNow and ShockNext are scalars: no further extraction needed
  
% JUMP VARIABLES:
  Vnow = xnow(1:nV);
  Vnext = xnext(1:nV);
  Cnow = xnow(nV+1);
  Cnext = xnext(nV+1);
  pnow = xnow(nV+2);
  pnext = xnext(nV+2);

% ENDOGENOUS STATE VARIABLES:
  PhiHatNow = xnow(nV+3:nV+2+nPhi);
  PhiHatNext = xnext(nV+3:nV+2+nPhi);
  
%LAST ELEMENT HAS BEEN LEFT IN, NO NEED TO REATTACH IT:
%% Attaching last element
%  PhiHatNow = [PhiHatNow; max(0,1-sum(PhiHatNow))];
%  PhiHatNext = [PhiHatNext; max(0,1-sum(PhiHatNext))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESHAPE VECTORS TO MATRICES size (nump,nums)
  [nump,nums]=size(PMAT);
  siz=[nump,nums];
  PhiHatNow = reshape(PhiHatNow,siz);
  PhiHatNext = reshape(PhiHatNext,siz);
  Vnow = reshape(Vnow,siz);
  Vnext = reshape(Vnext,siz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE SOME PRELIMINARY VARIABLES

% Money growth
  munow  = mu*exp(znow);     %this is consistent with ShockNow ~ N(0,s^2); then take logs below.   
  munext = mu*exp(znext);    %%%%%%%%%%%%%%%%%%%%%%% OJO
  
% Calculate the transition matrix Rnow
  Rnow = sparse(nump,nump); 
  nowoffset = log(munow)/pstep;
  if nowoffset==0
      Rnow = eye(nump);
  else
      Rnow(1,1:ceil(nowoffset)) = 1;             
      remoffset = nowoffset - floor(nowoffset);
      startfirstdiag = [max([1 ; -floor(nowoffset)])  max([1 ; 1+ceil(nowoffset)])];         
      endfirstdiag = [min([nump-1 ; nump-ceil(nowoffset)])  min([nump ; nump+floor(nowoffset)])]; 
      startsecdiag = [max([2 ; 1-floor(nowoffset)])  max([1 ; 1+ceil(nowoffset)])];               
      endsecdiag = [min([nump ; nump-ceil(nowoffset)+1])  min([nump ; nump+floor(nowoffset)])];   
      Rnow(startfirstdiag(1):endfirstdiag(1),startfirstdiag(2):endfirstdiag(2)) = ...
           remoffset*speye(nump-ceil(abs(nowoffset)));
      Rnow(startsecdiag(1):endsecdiag(1),startsecdiag(2):endsecdiag(2)) = ...
      Rnow(startsecdiag(1):endsecdiag(1),startsecdiag(2):endsecdiag(2)) + ...
           (1-remoffset)*speye(nump-ceil(abs(nowoffset)));
      Rnow(nump,nump+floor(nowoffset)+1:nump) = 1; 
  end

%WARNING!!! this appears to break down if nextoffset=0 exactly!! 
% Calculate the transition matrix Rnext
  Rnext = sparse(nump,nump); 
  nextoffset = log(munext)/pstep    ;
  if nextoffset==0
      Rnext = eye(nump);
  else
      Rnext(1,1:ceil(nextoffset)) = 1;             
      remoffset = nextoffset - floor(nextoffset)   ;
      startfirstdiag = [max([1 ; -floor(nextoffset)])  max([1 ; 1+ceil(nextoffset)])];         
      endfirstdiag = [min([nump-1 ; nump-ceil(nextoffset)])  min([nump ; nump+floor(nextoffset)])]; 
      startsecdiag = [max([2 ; 1-floor(nextoffset)])  max([1 ; 1+ceil(nextoffset)])];               
      endsecdiag = [min([nump ; nump-ceil(nextoffset)+1])  min([nump ; nump+floor(nextoffset)])];   
      Rnext(startfirstdiag(1):endfirstdiag(1),startfirstdiag(2):endfirstdiag(2)) = ...
           remoffset*speye(nump-ceil(abs(nextoffset)));
      Rnext(startsecdiag(1):endsecdiag(1),startsecdiag(2):endsecdiag(2)) = ...
      Rnext(startsecdiag(1):endsecdiag(1),startsecdiag(2):endsecdiag(2)) + ...
           (1-remoffset)*speye(nump-ceil(abs(nextoffset)));
      Rnext(nump,nump+floor(nextoffset)+1:nump) = 1    ;
  end
  %check = sum(Rnext)     

% real wage 
  wnow  = chi*pnow*Cnow^gamma;
  wnext = chi*pnext*Cnext^gamma;
  
% Payoff today
  Unow = Cnow*(pnow^epsilon)*(PMAT.^(1-epsilon)-wnow*sMAT.*PMAT.^(-epsilon));   
 
% MAX and optimal prices with quadratic spline on V: NOW
  Mnow = NaN*ones(1,nums); pStarNow = NaN*ones(1,nums);     
  for col=1:nums
      [maxV, maxind] = max(Vnow(:,col));     %CAN PROBABLY BE DONE SIMULTANEOUSLY WITHOUT LOOPING
      if (any(maxind==1) || any(maxind==nump))
          disp('ERROR!! encountered corner solution for p in dyneqklein.m!')
          keyboard
      end
      localx = Pgrid(maxind-1:maxind+1);
      localx2 = localx.^2;
      localV = Vnow(maxind-1:maxind+1,col);
      XMAT = [ones(3,1) localx localx2];
      betaquad = (XMAT'*XMAT)\XMAT'*localV;
      pStarNow(col) = -0.5*betaquad(2)/betaquad(3);
      Mnow(col)  = [1 pStarNow(col) pStarNow(col)^2]*betaquad;
  end

% MAX and optimal prices with quadratic spline on V: NEXT
  Mnext = NaN*ones(1,nums); pStarNext = NaN*ones(1,nums);     
  for col=1:nums
      [maxV, maxind] = max(Vnext(:,col));
      if (any(maxind==1) || any(maxind==nump))
            disp('ERROR!! encountered corner solution for p in dyneqklein.m!')
            keyboard
      end
      localx = Pgrid(maxind-1:maxind+1);
      localx2 = localx.^2;
      localV = Vnext(maxind-1:maxind+1,col);
      XMAT = [ones(3,1) localx localx2];
      betaquad = (XMAT'*XMAT)\XMAT'*localV;
      pStarNext(col) = -0.5*betaquad(2)/betaquad(3);
      Mnext(col)  = [1 pStarNext(col) pStarNext(col)^2]*betaquad;
  end

  % Calculate adjustment values Dnow and Dnext
  Dnow  = ones(nump,1)*Mnow - Vnow;   
  Dnext = ones(nump,1)*Mnext - Vnext; 
  
% Calculate probabilities of adjustment Lambdanow and Lambdanext    

LambdaNow = adjustment(adjtype, Vnow, Dnow, wnow, ksi, alpha, lbar);
[LambdaNext ExpectGAINSnext] = adjustment(adjtype, Vnext, Dnext, wnext, ksi, alpha, lbar);
ExpectGAINSnext = ExpectGAINSnext - Vnext;

% Calculate PhiTildeNow from PhiHatNow = Phi_t-1
  PhiTildeNow = Rnow*PhiHatNow*TRANSMAT';
  
% Calculate matrix P which indicates the split of adjusting firms' 
% mass between the two grid points around the optimal price  
  PmatNow = zeros(nump,nums);
  for col=1:nums                           
    OPTindHI=find(Pgrid>pStarNow(col),1);
    OPTindLO=max(OPTindHI-1,1);
    PmatNow(OPTindHI,col)= (pStarNow(col)-Pgrid(OPTindLO))/pstep;
    PmatNow(OPTindLO,col)= (Pgrid(OPTindHI)-pStarNow(col))/pstep;    
  end                                                    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE RESIDUALS FOR JACOBIAN

% Value function residual
  VResid = -Vnow + Unow + beta * (pnow/pnext)*(Cnext/Cnow)^(-gamma) * ...
      Rnext' * (Vnext + ExpectGAINSnext) * TRANSMAT ;

% Euler residual
  eulerResid = Cnow^(-gamma) - nu*pnow - (beta/munext)*(pnow/pnext)*Cnext^(-gamma); 

% Price residual
  priceResid = pnow - (sum(sum(PhiHatNext.*(PMAT.^(1-epsilon)))))^(1/(1-epsilon));
  
% PhiHatNow = dist of prices at time of production NOW (after shocks and adjustment)
  PhiResid = -PhiHatNext + (1-LambdaNow).*PhiTildeNow + PmatNow.*(ones(nump,nump)*(LambdaNow.*PhiTildeNow)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONVERT BACK TO VECTORS
  VResid = VResid(:);
  PhiResid = PhiResid(:);
  %PhiResid(end) = [];    %DONT DELETE LAST ELEMENT
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEMBLE RESID VECTOR:
  Resid = [VResid; eulerResid; priceResid; PhiResid];
   
