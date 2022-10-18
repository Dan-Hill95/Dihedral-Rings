% Written by Dan Hill (2021) - University of Surrey, 
% adapted from codes by David Lloyd - University of Surrey and
% Daniele Avitabile - Vrije Universiteit Amsterdam;
% see Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016.
%
%% Inputs
%           x=[x[0],...,x[N]]           - Initial guess for the matching condition
%                                         For 6|m, a good initial guess will be of
%                                         the form x=y*[1,...,1], for some small y.
%                                         For 3~|m, 2|m, a good initial guess will be of
%                                         the form x=y*[-1,2,2,-1,2,2...], for some small y.
%                                         For 3|m, 2~|m, a good initial guess will be of
%                                         the form x=y*[-1,1,1,-1,1,1...], for some small y.
%           p=[mu, nu, kappa, m, N+1]   - Initial parameters of system; 
%                                         mu is the bifurcation parameter, 
%                                         nu is the quadratic coefficient,
%                                         kappa is the cubic coefficient,
%                                         m is the dihedral symmetry,
%                                         N+1 is the dimension of the ODE system
%           Dir                         - Direction of parameter continuation;
%                                         must be either 'pl' (plus), 'mn' (minus),
%                                         'sp' (small plus), or 'sm' (small minus).
%
%% Purpose
% Solves a N'th order Galerkin system for the stationary 2D 2-3 Swift-Hohenberg
% equation via finite-difference methods and continues solutions in mu-parameter
% space. For U(r,theta) such that
%
% U(r,theta) = u[0](r) + 2*sum_{i=1}^{N} u[i](r)*cos(m*i*theta), 
% then each u[i](r) satisfies 
%
% 0= F[i](u[i]):= -(1+d^2_r + 1/r*d_r - (2*m*i/r)^2)^2 u[i] - mu*u[i] + f[i](U).
%
% After solving an algebraic matching condition solutions are then 
% continued in mu-parameter space.
%
%% Outputs
%           branch- [Step, -1, mu, EucNorm(u), L2Norm(u[0]), max(|u[0]|), ..., max(|u[N]|)]
%
% All data is stored in a folder named as "D[m]_Patch_[Dir]" 
%
function branch = Cont_Patch(x,p,Dir)

close all, clc;

m = p(4);                   % Dihedral lattice
n = p(5);                   % ODE dimension
k = length(x);              % Matching problem dimension

% If the dimension k of the matching problem is larger than the dimension n of
% the ODE system, we raise the truncation order N (=n-1).
if k > n
n = k;
end
% Otherwise, the k-dimensional matching solution is embedded in the 
% n-dimensional ODE system.

% Setting up mesh parameters, collected as mesh_params, including
% finite-difference matrices
SetupDiffMats_Patch;

%% Initial data

a = MatchSoln(x,m,L,p(1));  % Solving the algebraic matching condition.

% The elements of a are the coefficients of our localised radial amplitudes,
% and a top-down profile of the solution is plotted.

pause;
close all                   % Closing the plot form MatchSoln.

% u0 is the initial guess for fsolve in the k-dimensional ODE system
u0 = zeros(k*N,1);          
% Defining u0 as in the plot in MatchSoln for initial guess
for i=0:k-1
    u0(1+(i)*N:(i+1)*N)= a(i+1)*(-1)^(m*i)*besselj(m*(i),r).*exp(-sqrt(p(1))*r);
end
    
myproblemHandle = @(u) Equation_Patch(u,p,mesh_params);
    
% Fsolve options
options = optimset('Jacobian','on','Display','iter','MaxIter',50,'TolFun',1e-7,'DerivativeCheck','off');

%Solving the "k"'th-order Galerkin truncation of the Swift-Hohenberg
%equation close to the initial guess
[u_out,fval,exitflag,output,jacobian] = fsolve(myproblemHandle,u0,options);

% Embedding the k-dimensional solution in the "n"'th-order Galerkin truncation
u1 = [u_out; zeros((n-k)*N,1)];

%% Plotting the surface of the found solution

hfig1=figure;
t = 0:0.01:2*pi;                % Angular variable
t = t';
[R,T]=meshgrid(r,t);
[U,T0]=meshgrid(u1,t);
UU(:,:,1)=U(:,1:N)/2;
  for i=1:n-1
      UU(:,:,i+1)=U(:,1+i*N:(i+1)*N);
      UU(:,:,i+1)=UU(:,:,i+1).*cos(m*i.*T);
  end
z=surf(R.*cos(T), R.*sin(T),sum(UU,3));
view(0,90);
z.FaceColor='flat';
z.FaceAlpha=1;
z.EdgeColor='none';
axis([-0.5*L*sqrt(2) 0.5*L*sqrt(2) -0.5*L*sqrt(2) 0.5*L*sqrt(2) -1 2]);
  
pause;
close(hfig1);

%% Define handle to right-hand side and time output function
prob     = @(u,p) Equation_Patch(u,p,mesh_params);
plotSol  = @(u,p,parent) PlotSurface_Patch(u,p,parent,mesh_params);
solMeas  = @(step,u,p) SolutionMeasures_Patch(step,u,p,mesh_params);

% Default code does not compute the stability of the solution, for computational speed
compSpec = [];%@(u,p) ComputeEigenvalues_Patch(u,p,prob);
plotSpec = [];%@(d,p,parentHandle) PlotSpectrum_SH(d,p,parentHandle);

%% Assign problem 
stepperPars.iContPar      = 1;
% Define the direction and step distance for continuation
if Dir == 'pl'
stepperPars.s0            = 0.01;
stepperPars.sMin          = 1e-9;
stepperPars.sMax          = 0.05;
elseif Dir == 'mn'
stepperPars.s0            = -0.01;
stepperPars.sMin          = 1e-9;
stepperPars.sMax          = 0.05;
elseif Dir== 'sp' 
stepperPars.s0            = 1e-4;
stepperPars.sMin          = 1e-9;
stepperPars.sMax          = 1e-3;
elseif Dir== 'sm' 
stepperPars.s0            = -1e-4;
stepperPars.sMin          = 1e-9;
stepperPars.sMax          = 1e-3;
else
    error('Final argument must be "pl" (forwards), "mn" (backwards), "sp" (small values forwards), or "sm" (small values backwards).');
end

stepperPars.pMin          = 0;
stepperPars.pMax          = 2;
stepperPars.maxSteps      = 2000;
stepperPars.nPrint        = 1;
stepperPars.nSaveSol      = 1;
stepperPars.finDiffEps    = 1e-7;
stepperPars.fsolveOptions = optimset('Display','off',...
                                     'DerivativeCheck','off',...
                                     'TolFun',1e-7,...
                                     'Jacobian','on',...
                                     'MaxIter',10);
stepperPars.optNonlinIter = 10;
stepperPars.dataFolder    = ['D',num2str(m),'_Patch_' Dir];
stepperPars.PlotSolution  = plotSol;
stepperPars.BranchVariables = solMeas;
% stepperPars.PlotBranchVariableId = [];%2;
stepperPars.ComputeEigenvalues = compSpec;
stepperPars.PlotSpectrum = plotSpec;      
stepperPars.PlotBranchVariableId = 1;


branch = SecantContinuation(prob,u1,p,stepperPars,'Branch');
end