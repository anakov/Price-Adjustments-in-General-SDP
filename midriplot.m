% Plots the histogram of price changes for the AC Nielsen dataset of Midrigan
load acnielsen
niel=data;

lboundM=-0.5;
hboundM=0.5; 
edges=[-inf linspace(lboundM,hboundM,24) inf];                         % boundaries of price change intervals

nielN = histc(niel,edges);
nielN = nielN(1:end-1);    % the last bin counts any values that match EDGES(end); see help histc
nielN = nielN./sum(nielN);

step = (hboundM-lboundM)/(length(nielN)-1);

figure
colormap([0.73 0.83 0.96])
bar(lboundM:step:hboundM,nielN,1,'EdgeColor','none')
title('AC Nielsen')
xlabel('Size of price changes')
ylabel('Density of price changes')
xlim([-0.5 0.5])


