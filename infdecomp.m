% Decomposes inflation response into intensive, extensive, and selection effects

%The steady state objects needed for this calculation are:
%   Pdisteroded: beginning-of-period distribution
%   lambda: the whole function lambda(p,a), nump x nums
%   Pdistans: matrix of desired LOG price adjustments (pstar-Pgrid), nump x nums 

format compact

%period in which money shock occurs - must be same as in DISTSIM!
if exist('shocktime','var')
   d_mu_shock = d_mu_path(shocktime);
else
   d_mu_shock = 1/scalefactor;
end

if d_mu_shock == 0, d_mu_shock = 1/scalefactor; end     % if there is no shock, take differences


% INFLATION IMPACT RESPONSE
dPI_dmu = (PI_identity_mean(time)- PI_SS) / d_mu_shock; 
%dPI_dmu = (PI_path_mean(time)- PI_SS) / d_mu_shock;
%dPI_dmu = (PI_path(time)- PI_SS) / d_mu_shock;

% KLENOW-KRYVTSOV DECOMPOSITION
% inflation = freqpchanges*EPchange
dfreqpchanges_dmu = (frac_adjusters_path(time) - freqpchanges) / d_mu_shock;
dAvPchange_dmu    = (AvPchange_path(time) - EPchange) / d_mu_shock;

KKintensive = dAvPchange_dmu * freqpchanges;
KKextensive = dfreqpchanges_dmu * EPchange ;
KK_decomp = KKintensive + KKextensive ;

% INFLATION = sum(sum(Pdistans.*(lambda - sum(sum(lambda .*Pdisteroded))).*Pdisteroded)) + ...
%    sum(sum(Pdistans.*Pdisteroded ))  *  sum(sum(lambda.*Pdisteroded )) ;
    
selection = sum(sum(Pdistans.*(lambda - sum(sum(lambda .*Pdisteroded))).*Pdisteroded));
intensive = sum(sum(Pdistans.*Pdisteroded )) ; 
extensive = sum(sum(lambda.*Pdisteroded ))   ;

d_Pchanges = (Pchanges_path(:,:,time) - Pdistans) ./ d_mu_shock;
d_PhiTilde = (reshape(PhiTilde_path(:,time),nump,nums) - Pdisteroded) ./ d_mu_shock;

d_intens = sum(sum(d_Pchanges.*Pdisteroded + Pdistans.*d_PhiTilde));  % element by element so derivative okay
d_extens = (frac_adjusters_path(time) - freqpchanges) / d_mu_shock;

intensive_margin = d_intens*extensive; 
extensive_margin = d_extens*intensive; 
selection_effect = KK_decomp - extensive_margin - intensive_margin ;

