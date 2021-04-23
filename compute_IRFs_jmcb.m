% Computes dynamic paths as a function of initial state and shocks
% In particular computes impulse-response functions
% For simplicity, 'time' is same as vector index (i.e. starts at time=1)

% IN ORDER TO BE ABLE TO CONSIDER VERY SMALL SHOCKS, COMPUTE DEVIATIONS ONLY

% EXTRACTING VARIABLES DIRECTLY FROM KLEIN REPRESENTATION

% Extract components from DYNAMICS OF STATE
  d_mu_path = STATEHISTORY(end,:);                                         % money growth history
  mu_path = mu + d_mu_path;                                                % money growth history

  d_PhiHat_path = STATEHISTORY(1:end-1,:);                                 % history of LAGGED distributions
  PhiHat_path = Pdist(:)*ones(1,TT) + d_PhiHat_path;                       % history of LAGGED distributions
                                                                           % by definition, PhiHat_t = Phi_t-1
% Extract components from DYNAMICS OF JUMPS:
  d_V_path = JumpHistory(1:nV,:);                                          % value function history
  V_path = V(:)*ones(1,TT) + d_V_path;                                     % value function history

  d_C_path = JumpHistory(nV+1,:);                                          % consumption history
  C_path = Cbar + d_C_path;                                                % consumption history

  d_p_path = JumpHistory(nV+2,:);                                          % real aggregate price history
  p_path = pbar + d_p_path;                                                % real aggregate price history

  m_path = (1./p_path)/Avgmoney - 1;                                       % real money history (in % dev from SS)

% CONSTRUCT OTHER VARIABLES NOT INCLUDED IN KLEIN REPRESENTATION
  
  w_path = chi*p_path.*C_path.^gamma;                                      % real wage history
  R_path = (1 - nu*C_path.^gamma.*p_path).^(-12)-1;                        % annual net nominal interest rate path 

% MAX with quadratic spline on V
  Ms_path = NaN*ones(TT,nums); pStars_path = NaN*ones(TT,nums);     
  for time=1:TT
    for col=1:nums
      Vt=reshape(V_path(:,time),nump,nums);
      [maxV, maxind] = max(Vt(:,col));       
      if (any(maxind==1) || any(maxind==nump))
         disp('ERROR!! encountered corner solution for p in compute_IRFs.m!')
         keyboard
      end
         localx = Pgrid(maxind-1:maxind+1);
            %LOG of real price
         localV = Vt(maxind-1:maxind+1,col);
            %LEVEL of value
         XMAT = [ones(3,1) localx localx.^2];
         betaquad = (XMAT'*XMAT)\XMAT'*localV;
         pStars_path(time,col) = -0.5*betaquad(2)/betaquad(3);
            %LOG of real price
         Ms_path(time,col)  = [1 pStars_path(time,col) pStars_path(time,col)^2]*betaquad;
            %LEVEL of (maximum) value   
    end
  end
  
  D_path = kron(Ms_path',ones(nump,1)) - V_path;                             % gain from adjustment history OK
  Lambda_path = adjustment(adjtype, V_path, D_path, ones(gridsize,1)*w_path, ksi, alpha, lbar);

% Now construct PhiTilde (beginning-of-t distributions) from PhiHat (lagged distributions)
  PhiTilde_path = NaN*zeros(gridsize,TT);         
  Rmatrix_path = NaN*zeros(nump*nump,TT); 
  offset_path = NaN*zeros(1,TT);
  
  for time=1:TT
    Rnow = sparse(nump,nump); 
    nowoffset = log(mu_path(time))/pstep;
    if nowoffset==0
        Rnow = eye(nump);
        Rmatrix_path(:,time)=Rnow(:);
        PhiTildeNow = Rnow*reshape(PhiHat_path(:,time),nump,nums)*TRANSMAT';
        PhiTilde_path(:,time)= PhiTildeNow(:);       % PdistEroded history
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
        Rmatrix_path(:,time)=Rnow(:);
        PhiTildeNow = Rnow*reshape(PhiHat_path(:,time),nump,nums)*TRANSMAT';
        PhiTilde_path(:,time)= PhiTildeNow(:);       % PdistEroded history
    end
    offset_path(time) = nowoffset;
  end
  d_offset_path = offset_path - log(mu)/pstep;

% Now that we have constructed PhiTilde, we no longer need PhiHat.
  Phi_path = PhiHat_path;
  Phi_path(:,1) = [];  % here we are losing 1 period for Phi and everything computed from Phi
% So now I have history of Phi_t = distribution at time of production in period t.

  Adjusters_path = Lambda_path.*PhiTilde_path;
  frac_adjusters_path = sum(Adjusters_path);

  AvPchange_path = NaN*ones(1,TT);
  Yi             = NaN*ones(nump,nums,TT);
  labor_path     = NaN*ones(1,TT);
  Pchanges_path = NaN*ones(nump,nums,TT);

  for time=1:TT
    Pchanges_path(:,:,time) = ones(nump,1)*pStars_path(time,:) - Pgrid*ones(1,nums);
    AvPchange_path(time) = sum(sum(reshape(Pchanges_path(:,:,time),nump,nums)...
                           .*reshape(Adjusters_path(:,time),nump,nums)))./frac_adjusters_path(time);
    Yi(:,:,time) = (p_path(time)./PMAT).^(epsilon).*C_path(time);                   
  end

% p and C measures
  p_path_DS=NaN*ones(1,TT-1);
  C_path_DS=NaN*ones(1,TT-1);
  ElogP_path=NaN*ones(1,TT-1);   
  p_variance=NaN*ones(1,TT-1);  %RAW VARIANCE (not DS)
  delta_wedge=NaN*ones(1,TT-1);
  delta_wedge_fp=NaN*ones(1,TT-1);
  p_dispersion_DS=NaN*ones(1,TT-1);
  p_path_mean=NaN*ones(1,TT-1);   %RAW MEAN (not DS)
  Y_path_mean=NaN*ones(1,TT-1);   %RAW MEAN (not DS)

for t=1:TT-1
  p_path_DS(t) =  sum(Phi_path(:,t).*(PMAT(:).^(1-epsilon))).^(1/(1-epsilon));
  p_path_mean(t) =  sum(Phi_path(:,t).*PMAT(:));
  
  ElogP_path(t) =  sum(Phi_path(:,t).*log(PMAT(:)));
  C_path_DS(t)=sum((PMAT(:)./w_path(t)).^(1-epsilon).*Phi_path(:,t)).^(1/(gamma*(epsilon-1)))*chi^(-1/gamma); 

  p_variance(t) =  sum((log(PMAT(:))-ElogP_path(t)).^2.*Phi_path(:,t));     % variance (Woodford's)
  delta_wedge(t) =  sum(Phi_path(:,t).*sMAT(:).*(PMAT(:)./p_path(t)).^(-epsilon)); % relevant dispersion 
                     %(only one sum needed bc all vectorized)
 % flex price dispersion: 
  delta_wedge_fp(t) =  sum(sum(ergodicDistFlexPrice.*sMAT.*((ones(nump,1)*optflexprice)/p_path(t)).^(-epsilon)));
                     %(double sum needed bc matrix format)
  p_dispersion_DS(t) =  sum(Phi_path(:,t).*(PMAT(:).^(-epsilon)))*(p_path(t))^(epsilon);  
  
  Y_path_mean(t) = sum(sum( squeeze(Yi(:,:,t)).*reshape(Phi_path(:,t),nump,nums) )); 
  labor_path(t) = sum(sum( squeeze(Yi(:,:,t)).*sMAT.*reshape(Phi_path(:,t),nump,nums) ));                   
end

  %INFLATION MEASURES USING SIMPLE PRICE AGGREGATION
  PI_identity_mean = frac_adjusters_path.*AvPchange_path;
  
  if INITCONDIT==0
      initPI = mu-1; % steady-state inflation
  else
      initPI = PI_identity_mean(:,1); % inflation from the inflation identity: pi=sum(x*l*Phi)
  end

  PI_path_mean = [initPI p_path_mean(2:end)./p_path_mean(1:end-1).*mu_path(2:end-1)-1];  
  
  if mu>muUShigh || idioshocks   == -1 % this is better behaved at high inflation but misses last observation
     ex_ante_real_interest_rate = R_path(1:end-2)-PI_path_mean(2:end)-(Rbar-PI_SS);
  else           % this has last observation but is less acurate with high inflation
     ex_ante_real_interest_rate = R_path(1:end-1)-PI_identity_mean(2:end)-(Rbar-PI_SS);
  end

  if INITCONDIT==0
    % the time 2 money shock is unexpected in period 1, so expected
    % inflation equals steady-state
      ex_ante_real_interest_rate = [0 ex_ante_real_interest_rate(2:end)]; 
  end
  
  
% inflation decomposition
intensive_margin_path = zeros(1,TT-1); 
extensive_margin_path = zeros(1,TT-1);  
selection_effect_path = zeros(1,TT-1); 
for time=1:TT
  infdecomp;
  intensive_margin_path(time) = intensive_margin; 
  extensive_margin_path(time) = extensive_margin;  
  selection_effect_path(time) = selection_effect; 
end

%IRFchecks;


% CONSTRUCT MATRIX HISTORIES:
%   V_MATPATH = NaN*zeros([nump,nums,TT-1]);
%   PhiTilde_MATPATH = V_MATPATH;
%   PhiHat_MATPATH = V_MATPATH;
%   Phi_MATPATH = V_MATPATH;
%   for ttemp = 1:TT-1
%       V_MATPATH(:,:,ttemp) = reshape(V_path(:,ttemp),nump,nums);
%       PhiTilde_MATPATH(:,:,ttemp) = reshape(PhiTilde_path(:,ttemp),nump,nums);
%       PhiHat_MATPATH(:,:,ttemp) = reshape(PhiHat_path(:,ttemp),nump,nums);
%       Phi_MATPATH(:,:,ttemp) = reshape(Phi_path(:,ttemp),nump,nums);
%   end

  