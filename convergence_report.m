% Reports convergence progress

function Y = convergence_report(adjtype,finegrid,MMiter,sizeofMMgrid,pDIFF,piter,VDIFF,Viter,PdistDIFF,Pdistiter) 
clc
disp(sprintf('Model     : %d',adjtype))
disp(sprintf('Fine grid : %d',finegrid))
if ~isempty(MMiter)
   disp(sprintf('MM iter   : %d of %d',[MMiter, sizeofMMgrid]))
end
disp(sprintf('\n'))  
disp(sprintf('p DIFF    : %0.3d',pDIFF))
disp(sprintf('p iter    : %d',piter))
disp(sprintf('\n'))  
disp(sprintf('V DIFF    : %0.3d',VDIFF))
disp(sprintf('V iter    : %d',Viter))
%disp(sprintf('PDIFF     : %d',PdistDIFF))
%disp(sprintf('PDistiter : %d',Pdistiter))
Y=[];