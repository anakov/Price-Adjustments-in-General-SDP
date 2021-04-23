% Computes Kolmogorov-Smirnov statistic for equality of two cdf's
load acnielsen

% fit Midrigan'ss data into our histogram bins 
lbound = -0.5;
hbound =  0.5; 
edges  = [-inf linspace(lbound,hbound,24) inf];  % same bin edges as in calcstats; produces 25 bins

pdfdata = histc(data,edges);
pdfdata = pdfdata(1:end-1);                      % last bin counts any values that match EDGES(end); see help histc
pdfdata = pdfdata./sum(pdfdata);

step = (hbound-lbound)/(length(pdfdata)-1);
PriceChangeGrid = lbound:step:hbound;

s = length(data);
modeldata = [];
data = [];

% generate artificial data from the histogram "prob"
for i=1:length(prob)
   data_for_bin_i = ones(round(s*prob(i)),1)*PriceChangeGrid(i);
   modeldata = [ modeldata; data_for_bin_i];
end

% generate artificial data from the histogram "pdfdata"
for i=1:length(pdfdata)
   data_for_bin_i = ones(round(s*pdfdata(i)),1)*PriceChangeGrid(i);
   data = [ data; data_for_bin_i];
end

 [KS_reject_same ,P,KSSTAT] = kstest2(modeldata,data)
