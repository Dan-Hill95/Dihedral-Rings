% Adapted by Dan J Hill 2021
% Original Copyright 2016 by Daniele Avitabile
% See https://www.maths.nottingham.ac.uk/personal/pmzda/
%
% If you use this code, please cite
% Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016

function plotHandle = PlotSpectrum_SH(d,p,parentHandle)

   if isempty(parentHandle)
     scrsz = get(0,'ScreenSize');
     plotHandle = figure('Position',[3*scrsz(3)/4 scrsz(4)/2 scrsz(3)/4 scrsz(4)/4]);
     parentHandle = plotHandle;
    else
      plotHandle = parentHandle;
   end

   figure(parentHandle);
   cla, hold on;

   reSpan = [-2 4];
   imSpan = [-2 2];
   plot(real(d),imag(d),'.','MarkerSize',10);
   plot(zeros(1,100),linspace(imSpan(1),imSpan(2)));
   xlim(reSpan); ylim(imSpan); 
   grid on;
   hold off;
   
   %% Save
   print -dtiff spectrum.tiff

end
