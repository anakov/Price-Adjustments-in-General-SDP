% Computes impulse-responses based on Klein's state-space solution
%    then calls compute_IRFs.m to calc IRFs of all other variables.
%
% Jumps are: V, C, p. 
% States are: Phihat, z.  
%      Phihat is a predetermined, endogenous state; z is an exogenous state with shocks.
% All state and jump variables represent deviations from the steady state,
%      so STATEDYNAMICS must have all eigenvalues <= 1.
%
% Klein solution is:
%          states_t+1 = STATEDYNAMICS*states_t + [0 ; shocks_t+1]
%             jumps_t = JUMPS*states_t.
%
%for simplicity, 'time' is same as index of row vectors: starts at time=1.

% SPECIFY ZERO SHOCKS for the distribution
%zeroshocks = zeros(length(PdistVEC)-1,1);  %representing distribution without last element
zeroshocks = zeros(length(PdistVEC),1);    %INCLUDING LAST ELEMENT IN DISTRIBUTION

%SET TIME 1 DISTRIBUTION (both dist and shock can be 0 or nonzero at initial time)
if INITCONDIT==0   %STEADY STATE: this is the baseline case
   disp('DISTSIM') 
   disp('STARTING FROM STEADY STATE DISTRIBUTION')
   shiftdist = Pdist;
elseif INITCONDIT==1    %STARTING FROM SHIFTED DISTRIBUTION in p direction (like $ shock)
    shiftsteps = 2   %number of steps in p direction
    disp('DISTSIM') 
    disp(sprintf('START FROM SHIFTED DISTRIBUTION (%d steps in p direction):',shiftsteps))
    disp(sprintf('  equivalent to %0.3g pct INCREASE in money supply.',100*pstep*shiftsteps))
       %Note: an INCREASE in money supply is a DECREASE in real prices
    shiftdist = Pdist;
    if shiftsteps>= 0
        shiftdist(1,:) = sum(shiftdist(1:1+shiftsteps,:));
        shiftdist(2:nump-shiftsteps,:) = shiftdist(shiftsteps+2:nump,:);
        shiftdist(nump-shiftsteps+1:nump,:) = 0;
    elseif shiftsteps < 0
        error('ERROR! must program this!')
    end
elseif INITCONDIT==2    %STARTING FROM SHIFTED DISTRIBUTION in A direction (like tech shock)
    shiftsteps = 2   %number of steps in A direction
    disp('DISTSIM') 
    disp(sprintf('START FROM SHIFTED DISTRIBUTION (%d steps in A direction):',shiftsteps))
    disp(sprintf('  equivalent to %0.3g pct INCREASE in productivity.',100*pstep*shiftsteps))
           %Note: an INCREASE in productivity is a DECREASE in marginal cost
    shiftdist = Pdist;
    if shiftsteps>= 0
        shiftdist(:,1) = sum(shiftdist(:,1:1+shiftsteps)')'; %#ok<UDIM>
        shiftdist(:,2:nums-shiftsteps) = shiftdist(:,shiftsteps+2:nums);
        shiftdist(:,nums-shiftsteps+1:nums) = 0;
    elseif shiftsteps < 0
        error('ERROR! must program this!')
    end
elseif INITCONDIT==3    %STARTING ALL FROM SAME PRICE
    disp('DISTSIM') 
    disp(sprintf('ALL FIRMS START FROM SAME PRICE (middle point of price grid)'))
    shiftdist = zeros(size(Pdist));
    if rem(nump,2)==1        %nump is odd
       midp = ceil(nump/2)
       shiftdist(midp,:) = sum(Pdist);
    else                     %nump is even
        midp = nump/2
        shiftdist(midp,:) = sump(Pdist/2);
        shiftdist(midp+1,:) = sump(Pdist/2);
    end
elseif INITCONDIT==4    %ALL START FROM SAME PRICE and PRODUCTIVITY (like switching on the idiosyncratic shocks!!)
    disp('DISTSIM') 
    disp(sprintf('SWITCHING ON THE IDIO SHOCKS: all firms start from middle price, middle prod'))
    % perhaps should start from a distribution over prices for same productivity
    shiftdist = zeros(size(Pdist));
    if (rem(nump,2)==1 && rem(nums,2)==1)       %nump and nums are odd
       midp = ceil(nump/2)
       mids = ceil(nums/2)
       shiftdist(midp,mids) = sum(sum(Pdist));
    else           
       error('ERROR! must program this!')
    end
elseif INITCONDIT==5    %ALL START FROM OPTIMAL FLEXIBLE PRICE
    disp('DISTSIM') 
    disp(sprintf('SWITCHING ON THE STICKINESS: all firms start from OPTIMAL FLEXIBLE price'))
    Pflex=log(wbar)+sgrid+log(epsilon/(epsilon-1)); %#ok<NASGU>
    Pdisteroded = RECURSEMAT*Pdist*TRANSMAT';              % COMBINING TWO STEPS
    Pmoves=sum(ones(size(Pdisteroded)).*Pdisteroded);      % marginal density (over price) of firms adjusting price
    shiftdist = zeros(size(Pdist));
    for col=1:nums                                         % loop building the new joint density (p,mc)
    OPTindHI=find(Pgrid>Pflex(col),1);
    OPTindLO=max(OPTindHI-1,1);
    shiftdist(OPTindHI,col)=shiftdist(OPTindHI,col)+ ...
        Pmoves(col)*(Pflex(col)-Pgrid(OPTindLO))/pstep;    
    shiftdist(OPTindLO,col)=shiftdist(OPTindLO,col)+ ...
        Pmoves(col)*(Pgrid(OPTindHI)-Pflex(col))/pstep;   
    end 
elseif INITCONDIT==6    %ALL START FROM OPTIMAL STICKY PRICE, as if they had just been forced to adjust (EURO)
    disp('DISTSIM') 
    disp(sprintf('ONE FORCED ADJUSTMENT: all firms start from OPTIMAL STICKY price'))
    %THIS CODE ASSUMES CORRECT Pdist and Pstar are in memory.
    shiftdist = zeros(size(Pdist));
    %TEXT ADAPTED FROM GE_Pdist_iter:
    Pdisteroded = RECURSEMAT*Pdist*TRANSMAT';              % COMBINING TWO STEPS
    Pmoves=sum(ones(size(Pdisteroded)).*Pdisteroded);      % marginal density (over price) of firms adjusting price
    for col=1:nums                                         % loop building the new joint density (p,mc)
    OPTindHI=find(Pgrid>pstar(col),1);
    OPTindLO=max(OPTindHI-1,1);
    shiftdist(OPTindHI,col)=shiftdist(OPTindHI,col)+ ...
        Pmoves(col)*(pstar(col)-Pgrid(OPTindLO))/pstep;    
    shiftdist(OPTindLO,col)=shiftdist(OPTindLO,col)+ ...
        Pmoves(col)*(Pgrid(OPTindHI)-pstar(col))/pstep;   
    end 
else
    error('ERROR!! must specify initial conditions for simulation!')
end

if abs(sum(sum(shiftdist))-1 > 1e-10) 
    error('ERROR!! probabilities dont sum to 1!!!')
end
initdist = shiftdist - Pdist;   %deviation from steady state
initdist = initdist(:);


%PUT TOGETHER ALL TIME 1 INITIAL CONDITIONS
initSTATE = [initdist ; time1moneyShock];   %TIME 1 STATE
initJUMPS = JUMPS*initSTATE;

STATEHISTORY = NaN*zeros(nPhi+1,TT);
JumpHistory = NaN*zeros(nV+2,TT);
STATEHISTORY(:,1) = initSTATE;           %time 1 state
JumpHistory(:,1) = initJUMPS;            %time 1 jumps

%initialize time 1 state for loop:
STATENOW = initSTATE;

%NOW SIMULATE DYNAMICS
for time = 2:TT
   STATENOW = STATEDYNAMICS*STATENOW + [zeroshocks; moneyshocks(time)];
   JumpHistory(:,time)  = JUMPS*STATENOW;
   STATEHISTORY(:,time) = STATENOW;
end




