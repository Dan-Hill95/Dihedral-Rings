% Adapted by Dan J Hill 2021
% Original Copyright 2016 by Daniele Avitabile
% See https://www.maths.nottingham.ac.uk/personal/pmzda/
%
% If you use this code, please cite
% Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016
function [V,LAMBDA] = ComputeEigenvalues_Patch(uu,p,parent)

  %% Compute linear operators
  
  [~,J] = parent(uu,p);

  %% Call direct eigenvalue solver
  [V,LAMBDA] = eigs(full(J),1,0.1);
end
