% Plots steady-state distribution, histogram, value function and related objects

set(0,'DefaultFigureWindowStyle','docked')  % docks all figures

midriplot
hold on
stairs((lbound:step:hbound)-step/2,prob,'k','LineWidth',2')
title('Actual and simulated distribution of price changes')
legend('AC Nielsen','Model')

figure
plot(sgrid,optflexprice-log(wbar)-markup,':g')  % since marginal cost is really sgrid + log(wbar)
hold on
plot(sgrid,pstar-log(wbar)-markup,'r')
title('Optimal price policy')
xlabel('Log (inverse) productivity')
ylabel('Log target relative price')
axis tight

figure
if nums>1
   mesh(sgrid,Pgrid-log(wbar)-markup,D/MedV)
   xlabel('Log (inverse) productivity')
   ylabel('Log relative price')
 else
   plot(Pgrid-log(wbar)-markup,D/MedV)
   xlabel('Log relative price')
end
title('Adjustment gain')
axis tight

figure
if nums>1
   map=colormap;
   colormap([1 1 1; map])
   mesh(sgrid,Pgrid-log(wbar)-markup,PdistSmooth)
   xlabel('Log (inverse) productivity')
   ylabel('Log relative price')
else
   plot(Pgrid-log(wbar)-markup,Pdist)
   xlabel('Log relative price')
end
title('Stationary density of firms')
axis tight
%zl=zlim;

figure
%map=colormap;
%colormap([1 1 1; map])
if nums>1
   mesh(sgrid,Pgrid-log(wbar)-markup,Pdisteroded)
   xlabel('Log (inverse) productivity')
   ylabel('Log relative price')
else
   plot(Pgrid-log(wbar)-markup,Pdisteroded)
   xlabel('Log relative price')
end
title('Density of firms after mc shock and inflation')
axis tight
zl=zlim;
%zlim(zl)

figure
%map=colormap;
%colormap([1 1 1; map])
if nums>1
   mesh(sgrid,Pgrid-log(wbar)-markup,(1-lambda).*Pdisteroded)
   xlabel('Log (inverse) productivity')
   ylabel('Log relative price')
else
   plot(Pgrid-log(wbar)-markup,(1-lambda).*Pdisteroded)
   xlabel('Log relative price')
end
title('Density of non-adjusting firms')
axis tight
zlim(zl)

figure
%map=colormap;
%colormap([1 1 1; map])
if nums>1
   mesh(sgrid,Pgrid-log(wbar)-markup,lambda.*Pdisteroded)
   xlabel('Log (inverse) productivity')
   ylabel('Log relative price')
else
   plot(Pgrid-log(wbar)-markup,lambda.*Pdisteroded)
   xlabel('Log relative price')
end
title('Density of adjusting firms')
axis tight

figure
plot(Dgrid2/MedV*100,Lvalues2)
title('Lambda as a function of L (relevant range)')
xlabel('Loss from inaction (in % of firm''s median value)')
ylabel('Probability of adjustment')
hold on
axis tight
xlim([0 1])


figure
bar(Dgrid2(1:end-1)/MedV*100,diff(Lvalues2))
title('Lambda pdf')
xlabel('Loss from inaction (in % of firm''s median value)')
ylabel('Probability of adjustment')
xlim([0 1])
   
if CalvoMenuMetric<1-10*eps^.5
   if nums>1
      figure
      mesh(sgrid,Pgrid-log(wbar)-markup,lambdaSmooth)
      xlabel('Log (inverse) productivity')
      ylabel('Log relative price')
   else
      plot(Pgrid-log(wbar)-markup,lambda)
      xlabel('Log relative price')
   end
   title('Probability of adjustment Lambda')
   axis tight
end

% figure
% bar(0:npstep:1,density_lambda,1)
% title('Density of adjustment probabilities')
% xlabel('Adjustment probability lambda')
% ylabel('Density of firms (after shock)')
% xlim([-npstep/2 1+npstep/2])

figure
bar(0:Dstep:1,density_D,1)
title('Realized losses and hazard function')
xlabel('Non-adjustment loss L as % of median V')
ylabel('Density of firms')
hold on
plot(Dgrid2/MedV*100,Lvalues2,'r')
xlim([-Dstep/2 1+Dstep/2])
ylim([0 0.5])

figure
bar(0:Dstep:1,density_Derod,1)
title('Potential losses and hazard function')
xlabel('Non-adjustment loss L as % of median V')
ylabel('Density of firms')
hold on
plot(Dgrid2/MedV*100,Lvalues2,'r')
% plot(0:Dstep:1,Lvalues3.*density_Derod,'r')
xlim([-Dstep/2 1+Dstep/2])
ylim([0 1])
 
figure
if nums>1
   mesh(sgrid,Pgrid-log(wbar)-markup,V)
   xlabel('Log (inverse) productivity')
   ylabel('Log relative price')
 else
   plot(Pgrid-log(wbar)-markup,V)
   xlabel('Log relative price')
end
title('Value function')
axis tight

