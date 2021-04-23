% Computes Jacobian by forward differences
% Michael Reiter, Universitat Pompeu Fabra, September 2006
% Feel free to use, copy and modify at your own risk;
% Input:
% func: string, name of function
% x: point at which to take Jacobian
%    if function value at x is already known, then give x={starting point, function value}
% step: scalar, relative stepwidth;
% more arguments will be passed on to the function;
%
function jac = jacob_reiter(func,x,step,varargin)
  f0 = feval(func,x,varargin{:});
  set(0,'DefaultFigureWindowStyle','normal')       % undocks all figures
  disp(sprintf('Maximum residual: %d', max(abs(f0))))
  if max(abs(f0))>eps^.5
     disp('WARNING: LARGE STEADY-STATE ERROR!')
  end
  n = size(x,1);  
  m = size(f0,1); 
  disp(sprintf('\n'))  
  disp(sprintf('SIZE OF JACOBIAN: %0.1f million elements',m*n/1e6))  
  disp(sprintf('No. equations: %d,  No. arguments: %d',[m, n]))  
  jac = sparse([],[],[],m,n,round(0.05*m*n)); 
  x0 = x;
  h = waitbar(0,'Computing Jacobian. Please, wait...');
  for i=1:n
    if rem(i,100)==0, waitbar(i/n,h), end
%    step2 = step*max(1,abs(x0(i)));
    step2 = step;   %antons preferred alternative (there is a big difference for the menu cost model)
    x = x0;
    x(i) = x0(i) + step2;
    jac(1:m,i) = (feval(func,x,varargin{:}) - f0)/step2;
  end
  close(h) 
