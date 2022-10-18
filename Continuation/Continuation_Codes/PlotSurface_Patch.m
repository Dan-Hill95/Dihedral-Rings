% Copyright 2016 by Daniele Avitabile
% See https://www.maths.nottingham.ac.uk/personal/pmzda/
%
% If you use this code, please cite
% Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016

function plotHandle = PlotSurface_Patch(u,p,parentHandle,mesh_params)
r = mesh_params.r;
N = mesh_params.N;
m = mesh_params.m;

% n = mesh_params.n;
n = floor(length(u)/N);
   %% Position and eventually grab figure
   if isempty(parentHandle)
     scrsz = get(0,'ScreenSize');
     plotHandle = figure('Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/4 scrsz(4)/4]);
     parentHandle = plotHandle;
    else
      plotHandle = parentHandle;
   end

   figure(parentHandle);

   %% 1D plot
   
%    plot(r,u(1:N),'.-');
   
   %% 2D plot
   
   t = 0:0.01:2*pi;
   t = t';
  u = u(1:n*N);
  [R,T]=meshgrid(r,t);
  [U,T0]=meshgrid(u,t);
  UU(:,:,1)=U(:,1:N)/2;
  for i=1:n-1
      UU(:,:,i+1)= U(:,1+i*N:(i+1)*N);
      UU(:,:,i+1)= UU(:,:,i+1).*cos(m*i.*T);
  end
  z=surf(R.*cos(T), R.*sin(T),sum(UU,3));
  view(0,90);
  z.FaceColor='flat';
  z.FaceAlpha=1;
  z.EdgeColor='none';
   pbaspect([1 1 1]);
  axis([-50*sqrt(2) 50*sqrt(2) -50*sqrt(2) 50*sqrt(2) -1 2]);
   
   %% Save
   % print -depsc state.eps
   print -dtiff state.tiff

end

