% Computes the steady state distribution of firms on (price, productivity)
% version 17 april 2008, written for model detrended by M

PdistDIFF=inf;                                              % reset difference of Pdist
Pdistiter=0;                                                % reset counter

while (PdistDIFF>PdistDIFFtol && Pdistiter<1500)            % iterate to convergence of Pdist
    Pdistiter=Pdistiter+1;                                  % increment counter
    Pdisteroded = RECURSEMAT*Pdist*TRANSMAT';               % Distribution after MC shock and deflation
    Pmoves=sum(lambda.*Pdisteroded);                        % marginal density (over all prices) of adjusting firms 
    PdistNew=(1-lambda).*Pdisteroded;                       % joint density of firms not adjusting price
    for col=1:nums                                          % loop building the new joint density (p,mc)
  
  % Adjusting firms choose the optimal price which may lie inbetween nodes (obtained from the interpolation of V)  
  
    OPTindHI=find(Pgrid>pstar(col),1);                      % index of the grid point immediately above pstar
    OPTindLO=max(OPTindHI-1,1);                             % index of the grid point immediately below pstar
    PdistNew(OPTindHI,col)=PdistNew(OPTindHI,col)+ ...      % splitting of the mass between the points HI and LO:
        Pmoves(col)*(pstar(col)-Pgrid(OPTindLO))/pstep;     % assigning mass to point HI
    PdistNew(OPTindLO,col)=PdistNew(OPTindLO,col)+ ...      % assigning mass to point LO
        Pmoves(col)*(Pgrid(OPTindHI)-pstar(col))/pstep;   

    end                                                    
    
    PdistDIFF=gridsize*max(max(abs(PdistNew-Pdist)));      % sup norm normalized by gridsize
    Pdist=PdistNew;                                        % updating Pdist   
end