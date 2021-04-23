% Sets up matrices for general equilibrium model

%INITIALIZATION of matrices based on guess for pbar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Building Plongvec and pilongvec to nest the case of recursive policy 
if RECURSE==0
   Plongvec = Pgrid;   %in LOGS
elseif RECURSE==1
   Plongvec = kron(ones(numpi,1),Pgrid);   %REPEATS whole Pgrid vector vertically
   pilongvec = kron(pigrid,ones(nump,1));  %REPEATS each element of pigrid numP times
end

% Matrices representing idiosyncratic state: converted to LEVELS
sMAT=ones(numlong,1)*exp(sgrid);  % exogenous (INVERSE) idiosyncratic productivity with transition matrix TRANSMAT
PMAT=exp(Plongvec)*ones(1,nums);  % endogenous sticky idiosyncratic price

% Implied value of Cbar based on the guess for pbar
Cbar = ((1-beta/mu)/(nu*pbar))^(1/gamma);                 

% Implied value of the real PAYOFF based on the guess for pbar and Cbar
PAYOFFMAT = Cbar*(pbar^epsilon)*(PMAT.^(1-epsilon)-wbar*sMAT.*PMAT.^(-epsilon));

% Initial guess for value function     
V = PAYOFFMAT/(1-beta);                       % initialize it based on guess for PAYOFFMAT

% check for uniqueness with different intial values
Vlow  = 2*min(min(V))*ones(size(V));
Vhigh = 2*max(max(V))*ones(size(V));
% V = Vhigh; 


% Initial guess for the distribution of firms    
Pdist=zeros(numlong,nums);          
Pdist(ceil(numlong/2),ceil(nums/2))=1;   % initial probability distribution with unit atom in middle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%BUILDING TRANSMAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each ROW of TRANSMAT refers to some state at t+1, while each COL of TRANSMAT refers to some state at t.

sigma2_eps = stdMC^2*(1-rho^2);                 % variance of technology innovation  

if idioshocks==1
  [TRANSMAT, out1, out2, out3, out4, out5] = tauchen(rho,0,sigma2_eps,nums,numstddevs);

  TRANSMAT=TRANSMAT';                         % Productivity state transition matrix
  if any(abs(out1'-sgrid)>eps^.5), error('Grid mismatch'); end
  if any(abs(sum(TRANSMAT)-1)>eps^.5), error('A column of TRANSMAT does not sum to 1'); end

elseif idioshocks==0
  TRANSMAT=eye(nums);

elseif idioshocks==-1
  if rem(nums,2)==0, error('nums must be odd for rep agent case'), end    
  TRANSMAT=zeros(nums,nums);
  TRANSMAT((1+nums)/2,:)=ones(1,nums);
end

  
% COLUMNS of TRANSMAT should sum to one.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%STEADY STATE DISTRIBUTION OF PRODUCTIVITY SHOCKS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[shockvecs,shockeigs] = eig(TRANSMAT);
%[maxeig,maxindex] = max(diag(shockeigs));    %SHOULD BE 1
%SHOCKPROBSS = shockvecs(:,maxindex);         %EIGENVECTOR ASSOCIATED WITH EIGENVALUE 1
%checkshockss = [SHOCKPROBSS  TRANSMAT*SHOCKPROBSS]';
%this MAY be needed for simulation exercises. if not, delete it.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%BUILDING RECURSEMAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if RECURSE==0
   RECURSEMAT = sparse([],[],[],nump,nump,nump*2);                      % allocating memory 
   if mu==1
   RECURSEMAT=eye(nump);    
   % no recursive policy. so RECURSEMAT represents effects of inflation only
   else
   RECURSEMAT(1,1:ceil(offset)) = 1;                                   % empty if ceil(offset) <= 0
   remoffset = offset - floor(offset);
   startfirstdiag = [max([1 ; -floor(offset)])  max([1 ; 1+ceil(offset)])];              % rowstart, colstart
   endfirstdiag = [min([nump-1 ; nump-ceil(offset)])  min([nump ; nump+floor(offset)])]; % rowend, colend
   startsecdiag = [max([2 ; 1-floor(offset)])  max([1 ; 1+ceil(offset)])];               % rowstart, colstart
   endsecdiag = [min([nump ; nump-ceil(offset)+1])  min([nump ; nump+floor(offset)])];   % rowend, colend
   RECURSEMAT(startfirstdiag(1):endfirstdiag(1),startfirstdiag(2):endfirstdiag(2)) = ...
              remoffset*speye(nump-ceil(abs(offset)));
   RECURSEMAT(startsecdiag(1):endsecdiag(1),startsecdiag(2):endsecdiag(2)) = ...
              RECURSEMAT(startsecdiag(1):endsecdiag(1),startsecdiag(2):endsecdiag(2))+...
              (1-remoffset)*speye(nump-ceil(abs(offset)));
   RECURSEMAT(nump,nump+floor(offset)+1:nump) = 1;                    % empty if floor(offset) >= 0 
   end
%  full(RECURSEMAT)                                                      % just to check
elseif RECURSE==1 
   % RECURSIVE POLICY
   RECURSEMAT = sparse([],[],[],numlong,numlong,numlong*2);
   tempRMAT = sparse([],[],[],nump,nump,nump*2);
   for picount = 1:length(pigrid)
      tempRMAT = 0*tempRMAT;     %resetting
      infloffset = pigrid(picount)/pstep;         
      totoffset = offset - infloffset; 
      if totoffset==0 
         totoffset=round(totoffset);
         tempRMAT(1,1:totoffset) = 1;                                             % empty if offset negative
         startdiag = [max([1 ; 1-totoffset])  max([1 ; 1+totoffset])];            % rowstart, colstart
         enddiag = [min([nump ; nump-totoffset])  min([nump ; nump+totoffset])];  % rowend, colend
         tempRMAT(startdiag(1):enddiag(1),startdiag(2):enddiag(2)) = speye(nump-abs(totoffset));
         tempRMAT(nump,nump+totoffset+1:nump) = 1;
      else
         tempRMAT(1,1:ceil(totoffset)) = 1;                                            %empty if ceil(totoffset)<=0
         remoffset = totoffset - floor(totoffset);
         startfirstdiag = [max([1 ; -floor(totoffset)])  max([1 ; 1+ceil(totoffset)])];%rowstart, colstart
         endfirstdiag = [min([nump-1 ; nump-ceil(totoffset)])...
             min([nump ; nump+floor(totoffset)])];                                     %rowend, colend
         startsecdiag = [max([2 ; 1-floor(totoffset)])  max([1 ; 1+ceil(totoffset)])]; %rowstart, colstart
         endsecdiag = [min([nump ; nump-ceil(totoffset)+1])  min([nump ; nump+floor(totoffset)])]; %rowend, colend
         tempRMAT(startfirstdiag(1):endfirstdiag(1),startfirstdiag(2):endfirstdiag(2)) = ...
             remoffset*speye(nump-ceil(abs(totoffset)));
         tempRMAT(startsecdiag(1):endsecdiag(1),startsecdiag(2):endsecdiag(2)) = ...
            tempRMAT(startsecdiag(1):endsecdiag(1),startsecdiag(2):endsecdiag(2)) + ...
            (1-remoffset)*speye(nump-ceil(abs(totoffset)));
         tempRMAT(nump,nump+floor(totoffset)+1:nump) = 1;                          %empty if floor(totoffset) >= 0 
      end
      RECURSEMAT((picount-1)*nump+1:picount*nump,(picount-1)*nump+1:picount*nump) = tempRMAT;
   end
end
% RECURSEMAT is defined so that its COLUMNS sum to one
% Each ROW of RECURSEMAT refers to some price or policy at t+1, and 
% each COL of RECURSEMAT refers to some price or policy at t.   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('MMiter','var')
MMiter=[];  
sizeofMMgrid=[];
end


