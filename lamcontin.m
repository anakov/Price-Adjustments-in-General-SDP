% Interpolates step adjustment function for the menu cost model

function LAMBDA = lamcontin(D,alpha)

[nr,nc] = size(D);

%D should be NEAR ZERO along its principal diagonal 
%   and INCREASING as we move away from diagonal.
%Therefore should be CONVEX as we move along the p direction:
% D2D = diff(diff(D));
% if any(any(D2D<0))
%     error('ERROR!! D is not convex!!')
% end


%LINEAR INTERPOLATION OF D at 'midpoints' of grid:
Dinterp = interp1((1:nr)',D,(0.5:nr+0.5)','linear','extrap');


%NOW CONSIDER THE 'HALF-INTERVALS' on left and right of each grid point.
%DEFINE LAMBDA ON EACH 'HALF-INTERVAL' AS FRACTION OF THAT HALF INTERVAL
%   ON WHICH LAMBDA=1.


%ALLOCATE LAMBDA TO LEFT HALF-INTERVALS:
Dlb = Dinterp(1:nr,:);
Dub = D;
LAMlb = Dlb>alpha;   %lambda at left endpoint of half-interval
LAMub = Dub>alpha;   %lambda at right endpoint of half-interval

LAMleft = NaN*D;
LAMleft(LAMlb==1 & LAMub==1) = 1;
LAMleft(LAMlb==0 & LAMub==0) = 0;
CHANGES = LAMlb~=LAMub;

LAMleft(CHANGES) = LAMlb(CHANGES).*((alpha-Dlb(CHANGES))./(Dub(CHANGES)-Dlb(CHANGES))) + ...
                   LAMub(CHANGES).*((Dub(CHANGES)-alpha)./(Dub(CHANGES)-Dlb(CHANGES)));
               


%ALLOCATE LAMBDA TO RIGHT HALF-INTERVALS:
Dlb = D;
Dub = Dinterp(2:nr+1,:);
LAMlb = Dlb>alpha;   %lambda at left endpoint of half-interval
LAMub = Dub>alpha;   %lambda at right endpoint of half-interval

LAMright = NaN*D;
LAMright(LAMlb==1 & LAMub==1) = 1;
LAMright(LAMlb==0 & LAMub==0) = 0;
CHANGES = LAMlb~=LAMub;

LAMright(CHANGES) = LAMlb(CHANGES).*((alpha-Dlb(CHANGES))./(Dub(CHANGES)-Dlb(CHANGES))) + ...
                   LAMub(CHANGES).*((Dub(CHANGES)-alpha)./(Dub(CHANGES)-Dlb(CHANGES)));
               
               
%NOW ALLOCATE LAMBDA TO ORIGINAL GRID POINTS:
LAMBDA = 0.5*(LAMleft+LAMright);
%lambda at grid point is average of lambdas on left and right half-intervals.
