% Implements Klein's QZ decomposition method for solving linear RE models
% Based on Kleins gauss program solve.m, 
%    rewritten for MATLAB 3 March 2008 by Jim Costain.
%
% Note: this is an older version which involves solving a Sylvester equation 
% The latest version used in "Distributional dynamics" avoids this by rearranging 
% matrices A and B and is a bit faster (~5 seconds)
%
% Idea is to generalize Blanchard/Kahn to cases where their method fails
%    due to invertibility problems. Klein uses QZ decomposition instead
%    of eigendecomposition because it ensures invertibility while decomposing
%    into stable and unstable directions as in blanchard/kahn.
%
% Model analyzed is a*E_t(x_t+1) + b*x_t + c*E_t(z_t+1) + d*z_t = 0
%
%    where z_t+1 = phi*z_t + eps_t+1, and eps_t+1 is white noise.
%
%    Here x_t = [d_t ; k_t], where k_t is predetermined: k_t+1 = E_t(k_t+1).
%    User must specify nk, the number of predetermined variables,
%    and must order the x vector accordingly.
%
% Define w_t = [k_t ; z_t].  Klein finds the unique stable solution,
%    if it exists, as follows:
%
%    d_t = F*w_t    w_t+1 = P*w_t + [0 ; eps_t+1].
%
% Also write entire vector as v_t = [d_t ; k_t ; z_t].

function [F,P,stableeigs,unstableeigs] = kleinsolve(a,b,c,d,phi,nk,eqcutoff)

nx = size(a,1);
nz = size(phi,1);
nd = nx-nk;
nw = nk+nz;
nv = nx+nz;

A = [a    c ; zeros(nz,nx)  eye(nz)];
B = [-b  -d ; zeros(nz,nx)  phi    ];

save abcd a b c d
clear a b c d
% New representation is:  A*E_t[d_t+1 ; w_t+1] = B*[d_t ; w_t]

toc
disp('Matrices loaded. Start QZ decomposition.') % a.k.a. generalized Schur decomposition

A=full(A);
B=full(B);

[S,T,Q,Z] = qz(A,B);
save AB A B
clear A B
% This produces the QZ decomposition:
%   Q*A*Z = S --> Q*A = S*Z'
%   Q*B*Z = T --> Q*B = T*Z'
%   where S and T are upper triangular (or block upper triangular)
%   and Q and Z are unitary, which means Q*Q' = I, Z*Z' = I
%   where prime denotes transpose (if real) or conjugate transpose (if complex).

disp(sprintf('\n'))  
toc
disp('QZ decomposition finished. Start Ordeig.')

E = ordeig(T,S);    %Note T,S to compute T_ii/S_ii.
%E = eig(T,S);      %Note T,S to compute T_ii/S_ii. % On my laptop ORDEIG is not available 

disp(sprintf('\n'))  
toc
disp('Ordeig finished. Start OrdQZ.')

selection = abs(E)>1+eqcutoff;
clear E
[S,T,Q,Z] = ordqz(S,T,Q,Z,selection);
% This reorders matrices so that |T_ii/S_ii| > 1+eqcutoff at upper left, and
%                                |T_ii/S_ii| <= 1+eqcutoff at lower right. 

eigenvalues = diag(T)./diag(S);
%  Ratios |T_ii/S_ii| represents generalized eigenvalues of problem:
%     lambda*A*v = B*v  -->  lambda*Q*A_v = Q*B*v  -->  lambda*S*Z'*v = T*Z'*v
%  Define transformed variables y = z'*v
%     then generalized eigenvalue problem is lambda*S*y = T*y.
%  Note eivals lambda represent rate of increase of system y_t+1 = inv(S)*T*y_t.

unstableeigs = eigenvalues(abs(eigenvalues)>1+eqcutoff);
check_unstable = [nd  length(unstableeigs)];

stableeigs = eigenvalues(abs(eigenvalues)<=1+eqcutoff);
check_stable  = [nw  length(stableeigs)];

disp(sprintf('\n'))  
disp(sprintf('Number of Jumps : %d  Number of Unstable: %d',check_unstable))
disp(sprintf('Number of States: %d  Number of Stable  : %d',check_stable))

% NUMBER OF UNSTABLE EIGENVALUES SHOULD EQUAL nd:
if (length(unstableeigs)<nd || length(stableeigs)<nw)
    disp(sprintf('\n'))  
    disp('ERROR!!! Wrong eigenvalue count!!')
    keyboard
end

Suu = S(1:nd,1:nd);      
Tuu = T(1:nd,1:nd);      

Sus = S(1:nd,nd+1:nv);   
Tus = T(1:nd,nd+1:nv);   

Sss = S(nd+1:nv,nd+1:nv);
Tss = T(nd+1:nv,nd+1:nv);


% By definition, using fact that Z is unitary:
% [d;   =   [z11  z12; * [u;   -->   [u;   =   [z11'  z21'; * [d;
%  w]        z21  z22]    s]          s]        z12'  z22']    w]

% Unitary  -->  Z'*Z=I  -->  z11'*z12=-z21'*z22  -->  z12*inv(z22)=-inv(z11')*z21'
%   Also,       Z'*Z=I  -->  z12'*z12+z22'*z22=I

z11 = Z(1:nd,1:nd);
z12 = Z(1:nd,nd+1:nv);
z21 = Z(nd+1:nv,1:nd);
z22 = Z(nd+1:nv,nd+1:nv);

disp(sprintf('\n'))  
toc
disp('OrdQZ finished. Start z22 rank check ')
% Getting the rank of z22
[U_z22,SV_z22,V_z22] = svd(z22);   %SVD SHOULD GIVE U*SV*V' = z, with SV diag and U,V unitary.
svdprod1 = U_z22*SV_z22*V_z22';
CHECKSVD1 = max(max(abs(z22 - svdprod1)));

svdprod2 = z22*V_z22;
svdprod3 = U_z22*SV_z22;
CHECKSVD23 = max(max(abs(svdprod2 - svdprod3)));

%SVD is used to calculate rank. So checksvd must be small.
if CHECKSVD1>eps^.5 || CHECKSVD23>eps^.5
  disp(sprintf('\n'))  
  disp('Problem with check SVD:')
  disp(CHECKSVD1)
  disp(CHECKSVD23)
end

SV_z22 = diag(SV_z22);
tol_rank = max(size(z22)) * eps(max(SV_z22));
rank_z22 = sum(SV_z22 > tol_rank);

if rank_z22  < size(z22,1)
    disp(sprintf('\n'))  
    disp('FATAL ERROR!! Invertibility condition violated.')
    disp(sprintf('\n'))  
    keyboard
end

disp(sprintf('\n'))  
toc
disp('z22 rank check finished. Start Sylvester solution')

% NOW IMPOSE SADDLE PATH STABILITY USING THE QZ DECOMPOSITION.
% System can be rewritten as 
%       Q*A*E_t(v_t+1) = Q*B*v_t --> S*Z'*E_t(v_t+1) = T*Z'*v_t
%       --> E_t(y_t+1) = inv(S)*T*y_t,  where y_t = Z'*v_t
%
% Since S and T are triangular, eigenvalues of inv(S)*T are just T_ii/S_ii.
%
% Break into blocks: y_t = [u_t ; s_t]
%    where u_t corresponds to i such that |T_ii/S_ii| > 1 (unstable).
%    Likewise, using formula for inverse of upper triangular matrix:
%
% E_t [ u_t+1 ]  =  [inv(Suu)  -inv(Suu)*Sus*inv(Sss) ] * [Tuu  Tus] * [u_t]
%     [ s_t+1 ]     [  0                     inv(Sss) ]   [ 0   Tss]   [s_t]
%                =  [inv(Suu)*Tuu   inv(Suu)*(Tus-Sus*inv(Sss)*Tss)] * [u_t]
%                   [   0                            inv(Sss)*Tss  ]   [s_t]
%
% Note dynamics of s_t are decoupled from those of u_t, and involve stable
%   eigenvalues only (up to one, subject to numerical error), 
%   so lim j->infty E_t s_t+j is finite.
%
% Now suppose u depends linearly on s: u=Fs, so that u converges to a finite limit too.
%   Then we must have
%
% E_t u_t+1 = F E_t s_t+1 = F*inv(Sss)*Tss s_t = 
%     inv(Suu)*Tuu*F s_t + inv(Suu)*(Tus-Sus*inv(Sss)*Tss) s_t for ANY s_t.
%
% Therefore F must solve the Sylvester equation 
% F = inv(Tuu)*Suu*F*inv(Sss)*Tss - inv(Tuu)*(Tus-Sus*inv(Sss)*Tss).
%
% NOTE: THIS WHOLE LYAPUNOV STEP IS AVOIDED IN A LATER VERSION BY REARRANGING 
% OF MATRICES AS IN KLEIN'S ORIGINAL METHOD
%
% This mapping is stable, since inv(Sss)*Tss and inv(Tuu)*Suu both have
% stable eigenvalues only. Therefore it can be solved by iteration.
% Alternatively, it can be solved by dlyap.m from MATLAB control toolbox.

Tuui = eye(nd)/Tuu;
sdyn = Sss\Tss;

% Sylvester equation is solved by dlyap.m:
% F = AA*F*BB + CC;

AA   = Tuui*Suu;
BB   = sdyn;
CC   = -Tuui*(Tus-Sus*sdyn);

F_y = dlyap(AA,BB,CC); %SYLVESTER solution
clear AA BB CC

sw = z22+z21*F_y;
inv_sw = eye(nw)/sw;
compF = (z11*F_y + z12)*inv_sw ;
compP = sw*sdyn*inv_sw;

maximagF = max(max(abs(imag(compF))));
if maximagF>2*eps^.5
    disp(sprintf('\n'))  
    disp('WARNING! Large imaginary part observed in F:'), disp(maximagF)
end

maximagP = max(max(abs(imag(compP))));
if maximagP>2*eps^.5
    disp(sprintf('\n'))  
    disp('WARNING! Large imaginary part observed in P:'), disp(maximagP)
end

biggeststable = max(abs(eig(compP)));
if biggeststable > 1+eqcutoff
   disp(sprintf('\n'))  
   disp('ERROR!! large eigenvector classified as stable!!')
   keyboard
elseif biggeststable > 1
   % disp(sprintf('\n'))  
   % disp('Note: biggeststable infinitesimally greater than one: rescaling down.')
   compP=compP/biggeststable;
end

F = real(compF);
P = real(compP);

Pk = P(1:nk,:);
Pz = P(nk+1:end,:);

% KLEIN SETUP IS: a*E_t(x_t+1) + b*x_t + c*E_t(z_t+1) + d*z_t = 0
%   starting from arbitrary [k_t;z_t], which determines d_t, etc.
%
%  Here E_t x_t+1 = [d_t+1; = [F*P; * [k_t;    and x_t = [d_t; = [ F ; * [k_t; 
%                    k_t+1]     Pk]    z_t]               k_t]    I 0]    z_t]

% Similarly, A*E_t(v_t+1) = B*v_t from arbitrary [k_t;z_t],
%  but not from arbitrary v_t.

%IF SOLUTION IS CORRECT THE FOLLOWING MATRICES SHOULD BE IDENTICALLY ZERO:
load abcd 
CHECKKLEINSOLVE = a*[F*P ; Pk] + b*[F ; [eye(nk) zeros(nk,nz)]] + c*Pz + d*[zeros(nz,nk) eye(nz)];
clear a b c d 

load AB
CHECKKLEINSOLVE2 = A*[F*P ; P] - B*[F ; eye(nw)];
clear A B

CHECKKLEINSOLVE3 = S*Z'*[F*P ; P] - T*Z'*[F ; eye(nw)];

checkklein = max(max(abs(real(CHECKKLEINSOLVE))));
checkklein2 = max(max(abs(real(CHECKKLEINSOLVE2))));
checkklein3 = max(max(abs(real(CHECKKLEINSOLVE3))));

if (checkklein > 2*eps^.5 || checkklein2 > 2*eps^.5 || checkklein3 > 10*eps^.5 )
    disp(sprintf('\n'))  
    toc
    disp('POSSIBLE ERROR: Identities of Klein model not satisfied?')
    disp(checkklein);
    disp(checkklein2);
    disp(checkklein3);
else 
    disp(sprintf('\n'))  
    toc
    disp('Sylvester solution finished')
    disp(sprintf('\n'))  
    disp('Klein check: OK')
end

delete AB.mat
delete abcd.mat