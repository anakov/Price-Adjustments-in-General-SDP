function EKappa = expectmc(adjtype, D, wbar, ksi, alpha, lbar)
% Computes the expected menu cost in the stochastic menu cost model

L =  D./wbar;                                          % Lambda depends on D/w, value in units of labor time

[nump,nums] = size(D);
EKappa = NaN*ones(nump,nums);

minD = min(min(D));                                      
maxD = max(max(D));

Dgrid = linspace(minD,maxD,101);
KappaGrid = Dgrid(1:end-1)/wbar;

lambda= adjustment(adjtype, zeros(size(Dgrid)), Dgrid', wbar, ksi, alpha, lbar);
dlambda = diff(lambda);

for j=1:nums
    for i=1:nump
       lastKappa_index = find(KappaGrid < L(i,j),1,'last');  
       EKappa(i,j) = KappaGrid(1:lastKappa_index)*dlambda(1:lastKappa_index);
    end
end

    