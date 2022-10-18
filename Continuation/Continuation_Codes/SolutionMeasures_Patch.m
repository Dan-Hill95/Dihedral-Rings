% Copyright 2016 by Daniele Avitabile
% See https://www.maths.nottingham.ac.uk/personal/pmzda/
%
% If you use this code, please cite
% Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016

function F = SolutionMeasures_Patch(step,uu,pp,mesh_params)

  %% Rename parameters
     r = mesh_params.r;
     N = mesh_params.N;
    Lr = max(r);
    hr = abs(r(2)-r(1));
     n = floor(length(uu)/N);
     u = uu(1:n*N);
  %% Compute quadrature weights and l2norm
     w = ones(n*N,1); 
l2Norm = sqrt(sum( hr * w .* u.^2)/(2*Lr));
v = zeros(n,1);
for j=1:n
v(j) = max(abs(u(1+(j-1)*N:j*N)));
end

  %% Allocate
     F = zeros(1,n+1);

  %% Assign branch variables
     F = [l2Norm v'];
  
end
