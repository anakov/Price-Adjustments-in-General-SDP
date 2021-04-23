% Computes steady-state statistics 
Adjusters=lambda.*Pdisteroded;                                       % density of adjusting firms
Pdistans=ones(nump,1)*pstar-Pgrid*ones(1,nums);                      % matrix of distances from optimal price

AbsPdistans=abs(Pdistans);                                           % matrix of absolute price changes
Pincreases=Pdistans.*(Pdistans>eps^.5);                              % matrix of price increases
Pdecreases=Pdistans.*(Pdistans<-eps^.5);                             % matrix of price decreases

freqpchanges=sum(sum(Adjusters));                                    % percent of firms changing price in a period
massPriceIncreases=sum(sum(Adjusters.*(Pdistans>eps^.5)));           % mass of firms which increase their price   
massPriceDecreases=sum(sum(Adjusters.*(Pdistans<-eps^.5)));          % mass of firms which increase their price    

fracPriceDecr=1-massPriceIncreases/freqpchanges;                     % fraction of price decreases

EPchange=sum(sum(Pdistans.*Adjusters))/freqpchanges;                 % mean price change conditional on change
MeanAbsPchange=sum(sum(AbsPdistans.*Adjusters))/freqpchanges;        % mean abs price change conditional on change
EPincrease=sum(sum(Pincreases.*Adjusters))/massPriceIncreases;       % mean price increase
EPDecrease=sum(sum(Pdecreases.*Adjusters))/massPriceDecreases;       % mean price decrease
AvAbsDistFromPstar=sum(sum(abs(Pdistans).*Pdist));                   % average absolute distance from optimal price

STDpchange=sqrt(sum(sum((Pdistans-EPchange).^2.*Adjusters))...
                               /freqpchanges);                       % standard deviation of price changes

KurtosisPchange=(sum(sum((Pdistans-EPchange).^4.*Adjusters ))...
                               /freqpchanges)/STDpchange^4;          % kurtosis of price changes

inflation=(1+freqpchanges*EPchange)^3-1;                             % inflation identity
decompCheck=1200*( (1+inflation)^(1/3)-mu );                         % inflation decomposition check
wageCheck=(wbar-chi*pbar*Cbar^gamma)/wbar*100;                       % Check wage identity

AverD=sum(sum(D.*Pdist));                                            % mean loss of all firms 
StdD=sqrt(sum(sum((D-AverD).^2.*Pdist)));                            % standard deviation of the loss of all firms 

vecD=D(:);                                                           % vetorizes matrix 
vecV=V(:);                                                           % vetorizes matrix 
vecPdist=Pdist(:);                                                   % vetorizes matrix 
vecPdisteroded=Pdisteroded(:);                                       % vetorizes matrix
vecAdjusters=Adjusters(:);                                           % vetorizes matrix
vecAbsPdistans=AbsPdistans(:);                                       % vetorizes matrix 
vecPincreases=Pincreases(:);                                         % vetorizes matrix 
vecPdecreases=Pdecreases(:);                                         % vetorizes matrix 

AbsPch_dens=[vecAbsPdistans vecPdist];                               % matrix with AbsPdistans and their density
AbsPch_dens=sortrows(AbsPch_dens,1);                                 % sort by AbsPdistans ascending 
MedAbsPdistans=AbsPch_dens(find(cumsum(AbsPch_dens(:,2))>=.5,1),1);  % median absolute distance from optimal price

D_density=[vecD vecPdist];                                           % matrix with D's and their density
D_density=sortrows(D_density,1);                                     % sort by D's ascending 
MedD=D_density(find(cumsum(D_density(:,2))>=.5,1),1);                % median D

V_density=[vecV vecPdist];                                           % matrix with V and their density
V_density=sortrows(V_density,1);                                     % sort by V ascending 
MedV=V_density(find(cumsum(V_density(:,2))>=.5,1),1);                % median V

Pinc_density=[vecPincreases  vecAdjusters/...
              massPriceIncreases.*(vecPincreases>eps^5)];            % matrix of price increases and their density
Pinc_density=sortrows(Pinc_density,1);                               % sort by price increases ascending 
MedPincrease=Pinc_density(find(cumsum(Pinc_density(:,2))>=.5,1),1);  % median  price increase

Pdec_density=[vecPdecreases vecAdjusters/...
              massPriceDecreases.*(vecPdecreases<-eps^5)];           % matrix of price decreases and their density
Pdec_density=sortrows(Pdec_density,1);                               % sort by price decreases ascending 
MedPdecrease=Pdec_density(find(cumsum(Pdec_density(:,2))>=.5,1),1);  % median price decrease

EPrice=sum((PMAT(:)).*Pdist(:));                                     % mean price (levels) (NOT D-S AGGREGATE !!)

ElogPrice=sum(log(PMAT(:)).*Pdist(:));                               % mean log-level price (NOT D-S AGGREGATE !!)

VarPrice=sum((log(PMAT(:))-ElogPrice).^2.*Pdist(:));                 % std of prices (levels)

PriceDispersion=sum(sum( (PMAT/pbar).^(-epsilon).*Pdist ));          % Dixit-Stiglitz price dispersion measure
WeightPriceDispers=sum(sum( sMAT.*((PMAT/pbar).^(-epsilon)).*Pdist )); % Productivity-weighted D-S dispersion      

Y_mean=Cbar*PriceDispersion;                                         % Simple average consumption

NewPriceDeviations=ones(nump,1)*pstar-ones(nump,nums)*log(EPrice);   % deviations of new prices from mean
NewPriceIncreases=NewPriceDeviations.*(NewPriceDeviations > eps^.5); % deviations of new price increases from mean
STDpincrease=sqrt(sum(sum((NewPriceIncreases).^2.*...         
    Adjusters.*(Pdistans>eps^.5)))/massPriceIncreases);              % standard deviation of price increases

absP_density=[vecAbsPdistans vecAdjusters/freqpchanges];             % absolute price changes and their density
absP_density=sortrows(absP_density,1);                               % sort by price changes ascending 
MedAbsPchange=absP_density(find(cumsum(absP_density(:,2))>=.5,1),1); % median absolute price change

lbound = -0.5;
hbound = 0.5;

edges=[-inf linspace(lbound,hbound,24) inf];                         % boundaries of price change intervals
nbins=length(edges)-1;                                               % number of price bins
prob=zeros(nbins,1);                                                 % initialize mass in each bin with zero
for i=1:nbins                                                        % density in nbins price change intervals
    prob(i)=sum(vecAdjusters((Pdistans>edges(i)...
                            & Pdistans<=edges(i+1))))/freqpchanges;  % count of firms in price change interval k 
end

fracSmallChanges=sum(vecAdjusters((abs(Pdistans)<=0.05)))...
                              /freqpchanges;                         % fraction of small price changes

Rbar   = (1 - nu*Cbar.^gamma.*pbar).^(-12)-1;                        % annual net nominal interest rate
Rss    = mu/beta;                                                    % monthly gross nominal interest rate
PI_SS = mu-1;                                                        % inflation

Avgmoney = 1/pbar;                                                   % Average money holdings

AvRevenue = sum(sum(Cbar*(pbar^epsilon)*PMAT.^(1-epsilon).*Pdist));  % Average (and total) firms' revenue

Nbar = (pbar^epsilon)*Cbar*sum(sum(Pdist.*sMAT.*PMAT.^(-epsilon)));  % aggregate hours

U = Cbar^(1-gamma)/(1-gamma) - chi*Nbar + nu*log(Avgmoney); 

if exist('alpha','var')
% Menu Cost flow is freqpchanges*alpha*wbar.                         %as share of: 
   McostinRev=100*freqpchanges*alpha*wbar/AvRevenue;                 %  the flow of revenue of all firms
   McostinWageBill=100*freqpchanges*alpha*wbar/(Nbar*wbar);          %  the total wage bill 
   McostinCbar=100*freqpchanges*alpha/Cbar;                          %  consumption (units?)
end

varlambda=sum(sum((lambda-freqpchanges).^2.*Pdisteroded));           % Cross-sectional variance of lambda   
CalvoMenuMetric=varlambda/(freqpchanges*(1-freqpchanges));           % Compute Calvo-Menu cost metric

top_p_firms    = sum(Pdisteroded(end,1:end));                        % firms in top price bin
bottom_p_firms = sum(Pdisteroded(1,1:end));                          % firms in bottom price bin
top_MC_firms   = sum(Pdisteroded(1:end,end));                        % firms in top cost bin
bottm_MC_firms = sum(Pdisteroded(1:end,1));                          % firms in bottom cost bin
firms_hitting_boundaries = top_p_firms + bottom_p_firms + ...
                           top_MC_firms + bottm_MC_firms ;           % firms hitting boundaries of grid

% hazard histogram (based on Pdisteroded!)
edges=linspace(-eps^.5,1,101);                                       % boundaries of hazard intervals
npbins=length(edges)-1;                                              % number of hazard bins
density_lambda=zeros(npbins,1);                                      % initialize mass in each hazard bin with zero
for i=1:npbins                                                       % build density in hazard intervals
    density_lambda(i)=sum(vecPdisteroded((lambda>edges(i)...
                                         & lambda<=edges(i+1))));    % count of firms in hazard interval k 
end
%density_lambda=[sum(vecPdisteroded(lambda<=eps^.5));density_lambda]; % adding mass at zero
npbins=length(density_lambda);                                        % number of p bins
npstep=(1/(npbins-1));
if abs(sum(density_lambda)-1)>eps^.5, 
   disp('Density L does not sum to one'), sum(density_lambda)
end

minD=min(min(D));                                      
maxD=max(max(D));
Dcutoff=0.01*MedV;                                                   % very few firms should forgo such huge gains

if maxD>minD && Dcutoff>minD
% prepare Lambda as a function of D graph 
   Dgrid=linspace(minD,maxD,201);
   Dgrid2=linspace(minD,Dcutoff,301);
   Lvalues = adjustment(adjtype, zeros(size(Dgrid)), Dgrid', wbar, ksi, alpha, lbar);
   Lvalues2 = adjustment(adjtype, zeros(size(Dgrid2)), Dgrid2', wbar, ksi, alpha, lbar);
end
  
% losses histogram (based on Pdist!)
edges=linspace(-eps^.5,Dcutoff+eps^.5,31);                          % boundaries of D intervals
nDbins=length(edges)-1;                                              % number of D bins
%lv3step=length(Lvalues2)/(nDbins-1);
%Lvalues3 = Lvalues2(round(linspace(1,length(Lvalues2),nDbins)))';

density_D=zeros(nDbins,1);                                           % initialize mass in each D bin with zero
for i=1:nDbins                                                       % build density in hazard intervals
    density_D(i)=sum(vecPdist((D>edges(i) & D<=edges(i+1))));        % count of firms in D interval k 
end
%density_D=[sum(vecPdist(D<=eps^.5));density_D];                      % adding mass at zero
if abs(sum(density_D)-sum(vecPdist(D<=Dcutoff)))>eps^.5, 
   disp('Density D does not sum'), sum(density_D),
   sum(vecPdist(D<Dcutoff))
end
nDbins=length(density_D);                                            % number of D bins
Dstep=(1/(nDbins-1));
        
% losses histogram (based on PdistEeroded!)
density_Derod=zeros(nDbins,1);                                       % initialize mass in each D bin with zero
for i=1:nDbins                                                       % build density in hazard intervals
    density_Derod(i)=sum(vecPdisteroded((D>edges(i) & D<=edges(i+1))));  % count of firms in D interval k 
end
if abs(sum(density_Derod)-sum(vecPdisteroded(D<=Dcutoff)))>eps^.5, 
   disp('Density Derod does not sum'), sum(density_Derod),
   sum(vecPdisteroded(D<=Dcutoff))
end
maxmaxV=max(max(V));                                                 % maximum attainable V

ergodicDistFlexPrice=zeros(nump,nums);
optflexprice=log(exp(sgrid)*epsilon/(epsilon-1)*wbar);               % log of optimal flexible price
for col=1:nums 
    OPTflexIndHI=find(Pgrid>optflexprice(col),1);                    % index of the grid point immediately above 
    OPTflexIndLO=max(OPTflexIndHI-1,1);                              % index of the grid point immediately below 
    ergodicDistFlexPrice(OPTflexIndHI,col)=...
       sum(Pdist(:,col))*(optflexprice(col)-Pgrid(OPTflexIndLO))/pstep; % assigning mass to point HI
    ergodicDistFlexPrice(OPTflexIndLO,col)=...
       sum(Pdist(:,col))*(Pgrid(OPTflexIndHI)-optflexprice(col))/pstep; % assigning mass to point LO
end

profits = Cbar*(pbar^epsilon)*(PMAT.^(1-epsilon)-wbar*sMAT.*PMAT.^(-epsilon));
averprof=sum(sum(profits.*Pdist));

optprof = Cbar*(pbar^epsilon)*((exp(optflexprice)).^(1-epsilon)-wbar*exp(sgrid).*(exp(optflexprice)).^(-epsilon));
averoptprof = optprof*(sum(ergodicDistFlexPrice))';

optrevenu = Cbar*(pbar^epsilon)*((exp(optflexprice)).^(1-epsilon));
avoptrev = optrevenu*(sum(ergodicDistFlexPrice))';

loss1 = (averoptprof - averprof)/ averoptprof;
loss2 = (averoptprof - averprof)/ avoptrev;

sail=(diag(lambda(OPTindex',1:nums)));
smoothsail=hpfilter(sail,1000);
lambdaSmooth=lambda;
for i=1:nums
    lambdaSmooth(OPTindex(i),i)=smoothsail(i);
end


sail=(diag(Pdist(OPTindex',1:nums)));
smoothsail=hpfilter(sail,1000);
PdistSmooth=Pdist;
for i=1:nums
    PdistSmooth(OPTindex(i),i)=smoothsail(i);
end
