% Solves dynamic general equilibrium

global Params;
format long

jacstep = 100*ptol; 
eqcutoff = 100*jacstep;

nV = nump*nums;
nPhi = nV;   

Ntot = nV+2+nPhi;    % variables in X vector are: V, C, p, Pdist
nz = 1;              % one exogenous driving process (money)

%NOW STORE VARIABLES IN PARAMS
Params.nV = nV; 
Params.nPhi = nPhi;
Params.nz = nz;

%Y WILL CONSIST OF FOUR PARTS
Params.ix = 1:Ntot;                   %length Ntot
Params.ixnext = Ntot+1:2*Ntot;        %length Ntot
Params.iz = 2*Ntot+1:2*Ntot+nz;             %length nz
Params.iznext = 2*Ntot+nz+1:2*Ntot+2*nz;    %length nz

Params.beta = beta; 
Params.gamma = gamma; 
Params.chi = chi; 
Params.epsilon= epsilon; 
Params.mu = mu; 
Params.nu = nu; 

Params.adjtype = adjtype; 
Params.lbar = lbar; 
Params.ksi = ksi;              
Params.alpha = alpha;          
Params.Calvoparam = Calvoparam; 

Params.PMAT = PMAT; 
Params.sMAT = sMAT; 
Params.TRANSMAT = TRANSMAT; 
Params.Pgrid = Pgrid; 
Params.pstep = pstep; 

PdistVEC = Pdist(:);

Xss = [V(:);Cbar;pbar;PdistVEC];
   %variables in X vector are: V, C, p, Pdist

stst2=[Xss;Xss;zeros(nz,1);zeros(size(nz,1))];

resid = feval(@dyneqklein_jmcb,stst2);
check = max(abs(resid));
   if(check>eps^.5 || 10*check>jacstep)
     disp('WARNING: Large residual at stst:');
     disp(check);
   end;
tic
disp(sprintf('\n'))  
disp('START JACOBIAN CALCULATION')


   jac = jacob_reiter(@dyneqklein_jmcb,stst2,jacstep);
   ajac = jac(:,Params.ixnext);
   bjac = jac(:,Params.ix);
   cjac = jac(:,Params.iznext);
   djac = jac(:,Params.iz);

disp(sprintf('\n'))  
disp('START KLEIN SOLUTION')

[JUMPS,STATEDYNAMICS,stableeigs,unstableeigs] = kleinsolve_jmcb(ajac,bjac,cjac,djac,phiz,nPhi,eqcutoff);

