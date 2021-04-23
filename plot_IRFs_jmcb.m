% Plots impulse-response functions
set(0,'DefaultFigureWindowStyle','docked')  % docks all figures

if adjtype==0          
   linescolor = 'b.-' ;
   MarkerSize = 5;
elseif adjtype==2          
   linescolor = 'rs-' ;
   MarkerSize = 3;
elseif adjtype==3          
   linescolor = 'gx-';
   MarkerSize = 6;
elseif adjtype==4          
   linescolor = 'mo-';
   MarkerSize = 5;
 elseif adjtype==5          
   linescolor = 'k^-';
   MarkerSize = 5;
end

figure(5)
subplot(1,2,1)
hold on
plot(scalefactor*(PI_path_mean-PI_SS),linescolor,'MarkerSize', MarkerSize)
title('Inflation')
subplot(1,2,2)
hold on
plot(scalefactor*(C_path/Cbar-1),linescolor,'MarkerSize', MarkerSize) % this is real consumption which enters utility
title('Consumption')
