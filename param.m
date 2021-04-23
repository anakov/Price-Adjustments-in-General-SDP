% Sets model parameters and sets up grid
% param.m: version 18 December 2008

% PROGRAM EXECUTION PARAMETERS                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idioshocks   = 1;                              % idiosyncratic shocks: 1-idio; 0-fixed heterogeneity; -1 rep agent
showconverge = 1*(version<3);                  % show convergence progress. set to 0 for calibrate
RECURSE      = 0;                              % 0 sticky prices only ; 1 sticky recursive plans 

% DISCOUNT AND MONEY GROWTH PARAMETERS         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r        = 1.04^(1/12)-1;                      % monthly net real interest rate
beta     = 1/(1+r);                            % utility discount factor
phiz     = 0;                                  % persistence of money growth shocks (Klein dynamic representation)
muACNiel = 0.999831;                           % AC Nielsen inflation
muDomin  = 1.002822;                           % Dominick's inflation
muGLlow  = 1.0064^(1/3);                       % GL07 baseline monthly gross money growth rate (2.58% annual)
muGLhigh = 1.21^(1/3);                         % GL07 high   monthly gross money growth rate (114.36% annual) 
muGanLow = 1.0037;                             % Gagnon low     (4.53% annual)
muGanMed = 1.0214;                             % Gagnon medium (28.93% annual)
muGanHi  = 1.0416;                             % Gagnon high   (63.08% annual)
muUShigh = 1.10^(1/12);                        % US high in the 1970s (10% annual)
mu       = 1;                                  % money growth used in this program run

% PREFERENCE PARAMETERS                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma   = 2;                                   % CRRA coefficient
epsilon = 7;                                   % price elasticity of demand
chi     = 6;                                   % labor supply parameter
nu      = 1;                                   % money demand parameter      
wbar    = chi*(1-beta/mu)/nu;                  % steady-state wage

% INITIAL GUESS FOR PRICE LEVEL                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
guess_pbar;

% estimate.m will search over:
% lbar: Calvo parameter
% alpha: fixed menu cost
% ksi: state dependence parameter
% rho: persistence of idiosyncratic productivity
% stdMC: standard deviation of log of idiosyncratic productivity
% NOTE: stdMC = (sigma2_eps/(1-rho^2))^.5

% SHOCK AND ADJUSTMENT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch adjtype
   case 0  % SSDP 0. SSDP 0. SSDP 0. SSDP 0. SSDP 0. SSDP 0. SSDP 0. SSDP 0. SSDP 0. SSDP 0. SSDP 0. SSDP 
     if version < 3
        if finegrid==0
        % best estimate on COARSE grid
        lbar  = 0.108928579454875;
        alpha = 0.031083222345823;
        ksi   = 0.293675662699625;
        rho   = 0.881156092163841;
        stdMC = 0.147421080785644;
       elseif finegrid==1
        % best estimate on FINE grid, distance = 0.057679501846595
        lbar  = 0.109100439345311;
        alpha = 0.031046390129951;
        ksi   = 0.290016668050705;
        rho   = 0.880791169770882;
        stdMC = 0.147714569025461;
        end
     end
   case 2 % CALVO 2. CALVO 2. CALVO 2. CALVO 2. CALVO 2. CALVO 2. CALVO  2. CALVO  2. CALVO  2. CALVO  2. CALVO 
        lbar  =  Calvoparam;
        alpha = [];
        ksi   = [];
        if version==1            
        % re-estimated on 25x25  
        rho   = 0.857625081201242;
        stdMC = 0.165363367309271;
        elseif version==2          
        % fixed productivity process
        rho   = 0.881156092163841;
        stdMC = 0.147421080785644;
        end
   case 3 % MENU COST 3. MENU COST 3. MENU COST 3. MENU COST 3. MENU COST 3. MENU COST 3. MENU COST 3. MENU COST
        if version==1                               
        % re-estimated on 25x25
        lbar  = [];
        ksi   = [];
        alpha =  0.0630548117601; % baseline
        rho   =  0.8468875723353;
        stdMC =  0.1445805243468;
        elseif version==2                        
        % fixed productivity process
        lbar  = [];
        alpha = 0.059080072006615;
        ksi   = [];
        rho   = 0.881156092163841;
        stdMC = 0.147421080785644;
        end
   case 4 % WOODFORD 4. WOODFORD 4. WOODFORD 4. WOODFORD 4. WOODFORD 4. WOODFORD 4. WOODFORD 4. WOODFORD 
        if version==1                              
        % re-estimated on 25x25,  distance = 0.078804403745638
        lbar  = 0.094587156333389;
        alpha = 0.062586215733648;
        ksi   = 1.334628921440485;
        rho   = 0.859568936509515;
        stdMC = 0.180462107484153;
        elseif version==2                         
        % fixed productivity process
        lbar  = 0.088711110455847;
        alpha = 0.000001824078809;
        ksi   = 1.748877395928954;
        rho   = 0.881156092163841;
        stdMC = 0.147421080785644;
        end
  case 5  % DKW 5. DKW 5. DKW 5. DKW 5. DKW 5. DKW 5. DKW 5. DKW 5. DKW 5. DKW 5. DKW 
     if version < 3
        if finegrid==0
        % best estimate on COARSE grid, distance = 0.050336599458295
        lbar  = 0.108579977283206;
        alpha = 0.034863241298866;
        ksi   = 0.231778803383148;
        rho   = 0.918308821866637;
        stdMC = 0.160007325336373;
        end
     end
   
end

% Build COST and PRICE GRIDS in LOGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if finegrid         
   nums = 101;                                   % number of cost grid points 
   NumPpoints = 101;                             % number of price points (approximate); nump will be the actual
   numstddevs=5;                                 % number of standard deviations around midpoint
else
   nums = 25;                                    % number of cost grid points 
   if accuracy==0
     NumPpoints = 25;                              % number of price points on coarse grid
   else
     NumPpoints = 31;                              % number of price points on coarse grid
   end
   numstddevs = 2.5;                             % number of standard deviations around midpoint
end
SMAX=numstddevs*stdMC;                           % maximum cost
SMIN=-numstddevs*stdMC;                          % minimum cost
sstep=(SMAX-SMIN)/(nums-1);                      % distance between grid points
sgrid=SMIN:sstep:SMAX;                           % marginal cost grid (log)     

markup = log(epsilon/(epsilon-1));                 % optimal flexible price markup

if idioshocks == -1,                               % representative agent case
   nums=1;                                         % single productivity state
   SMAX=0;                                         % single grid point at zero log productivity
   SMIN=0;                                         % productivity=1 in levels
   sgrid=0;
   sstep=0;  
   PMAX = SMAX + gridSpread;          
   PMIN = SMAX - gridSpread;          
%   PMAX = markup + log(wbar) + SMAX + gridSpread;          
%   PMIN = markup + log(wbar) + SMAX - gridSpread;          
else
%  Real price grid (LOGS) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   PMAX = markup + log(wbar) + SMAX ;              % maximum price (based on flex price policy)
   PMIN = markup + log(wbar) + SMIN ;              % minimum price
   extraspread = gridSpread*(PMAX-PMIN);           % stretching out grid for high inflation cases
   PMAX = markup + log(wbar) + SMAX + extraspread; % maximum price (based on flex price policy)
   PMIN = markup + log(wbar) + SMIN - extraspread; % minimum price
end

pstep = (PMAX-PMIN)/(NumPpoints-1);             % NumPpoints price intervals

offset = log(mu)/pstep;                         % noninteger number of price steps caused by inflation

Pgrid=(PMIN:pstep:PMAX)';                         % price grid (log)
nump=length(Pgrid);                               % number of price grid points

% Define INFLATION GRID in case of recursive policy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if RECURSE==0
   numlong = nump;
else
   %since Pgrid is defined in terms of log P, should define pigrid as CHANGES IN LOG P.
   pirad = 12;                                    % radius of inflation rates considered
   pistep = 0.002;                              % inflation step is half a pct point
   pigrid = log(mu)-pistep*pirad:pistep:log(mu)+pistep*pirad;
   numpi = length(pigrid);
   numlong = nump*numpi;
end

gridsize=numlong*nums;                            % total number of grid points 

% Convergence tolerance levels % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ctol = (10^(8-accuracy))*eps;                     % convergence tolerance for C iteration
ptol = (10^(8-accuracy))*eps;                     % convergence tolerance for p iteration         
Vtol = (10^(8-accuracy))*eps;                     % convergence tolerance for Value function iteration
PdistDIFFtol = (10^(8-accuracy))*eps;             % convergence tolerance for Pdist function iteration
warning('off', 'MATLAB:divideByZero')    
  