%Written by Dan Hill (2020) - University of Surrey

%Inputs: uu is all of our solutions, p is all of our parameters in the problem,
%mesh_params are the parameters of the mesh

%Coupled system of 4'th-order Galerkin truncation of the 2D
%Swift-Hohenberg equation in polar coordinates. For u(r,theta),

%0= F(u):= -(1+d^2_r + 1/r*d_r + 1/r^2*d^2_theta)^2 u - mu*u + a1*u^2 + a4*u^3.

% Solutions take the form u(r,theta) = u0(r) + 2*sum_{i=1}^{n-1} ui(r)*cos(2*m*i*theta)

%By projecting onto each Fourier mode, we arrrive at n equations to
%solve, found below.

%Outputs: F - Coupled Swift-Hohenberg system output for a given solution
%uu, J - Jacobian of the function F.
function [F,J] = Equation_Patch(uu,p,mesh_params)

%Calling parameters
N = mesh_params.N;          %length of the mesh
LN = mesh_params.LN;        %-(1+d^2_r + 1/r*d_r)^2 finite difference matrix
LM = mesh_params.LM;        %-(1+d^2_r + 1/r*d_r - (2*m*i)^2/r^2)^2 finite difference matrix for each u[i](r)
n = floor(length(uu)/N);    %Truncation order n
LM = LM(1:(n-1)*N,1:(n-1)*N);

u = uu(1:n*N);          %uu might hold extra information. We don't want that information when using u, hence we omit those components of uu

e0 = sparse(1:N,1:N,ones(1,N),n*N,n*N);


M1= sparse(n*N,n*N);
M2= sparse(n*N,n*N);
for j=0:n-1
M1 = M1 + sparse(1+j*N:n*N,1:(n-j)*N,repmat(u(1+j*N:(j+1)*N),n-j,1),n*N,n*N);
end
for j=0:n-1
M2 = M2 + sparse(1:(n-j)*N,1+j*N:n*N,repmat(u(1+(n-j-1)*N:(n-j)*N),n-j,1),n*N,n*N);    
end

%renaming parameters coming from p
mu = p(1);

%Defining the linear operator L:=-(1+d^2_r + 1/r*d_r - (2*m*j)^2/r^2)^2 - mu,
%for each j = 1, 2,..., n-1.
Lap = sparse(n*N,n*N);

LAP = blkdiag(LN, LM);
Lap = LAP - mu*speye(n*N);

%Pre-defining the size of variables, for speed
F = sparse(n*N,1);
Q = sparse(n*N,1);
C = sparse(n*N,1);

J0 = flip(speye(n*N));
for j=0:n-1
J0(1+j*N:(j+1)*N,1+(n-j-1)*N:(n-j)*N) = flip(J0(1+j*N:(j+1)*N,1+(n-j-1)*N:(n-j)*N));
end

% Quadratic nonlinearity

Q = 2*M2*J0*u + M1*u - 2*sparse(1:n*N,1:n*N,repmat(u(1:N),n,1),n*N,n*N)*u;

% Cubic nonlinearity

C = (M2*J0 + M1.' + M1 - sparse(1:n*N,1:n*N,repmat(u(1:N),n,1),n*N,n*N))*Q ...
        + ( (M2.')*M2 - sparse(1:n*N,1:n*N,repmat(Q(1+(n-1)*N:n*N),n,1),n*N,n*N)*J0 ...
            - sparse(1:n*N,1:n*N,repmat(Q(1:N),n,1),n*N,n*N))*u;

% F will be the final output function
F = Lap*u + p(2).*Q + p(3).*C;

if nargout > 1 %Calculating the Jacobian
  
J = sparse(n*N,n*N);
DQ = sparse(n*N,n*N);
DC = sparse(n*N,n*N);

Mu2=M2*u;

U = sparse(n*N,N);
MQ1= sparse(n*N,n*N);
M21= sparse(n*N,n*N);
MQ2= sparse(n*N,n*N);
for j=0:n-1
MQ1 = MQ1 + sparse(1+j*N:n*N,1:(n-j)*N,repmat(Q(1+j*N:(j+1)*N),n-j,1),n*N,n*N);
M21 = M21 + sparse(1+j*N:n*N,1:(n-j)*N,repmat(Mu2(1+j*N:(j+1)*N),n-j,1),n*N,n*N);
U  = U + sparse(1+j*N:(j+1)*N,1:N,u(1+j*N:(j+1)*N),n*N,N);
end
for j=0:n-1
MQ2 = MQ2 + sparse(1:(n-j)*N,1+j*N:n*N,repmat(Q(1+(n-j-1)*N:(n-j)*N),n-j,1),n*N,n*N);    
end

% Quadratic nonlinearity

DQ = 2*M2*J0 + 2*M1 - 2*sparse(1:n*N,1:n*N,repmat(u(1:N),n,1),n*N,n*N) ...
        + 2*M1.' - 2*M1*e0;

% Cubic nonlinearity

DC= M2*J0*DQ + MQ1.' + M1*DQ + MQ1 + M1.'*DQ + MQ2*J0 ...
    - sparse(1:n*N,1:n*N,repmat(u(1:N),n,1),n*N,n*N)*DQ - MQ1*e0 + 2*(M2.')*M2 ...
    + M21*J0 - sparse(1:n*N,1:n*N,repmat(Q(1+(n-1)*N:n*N),n,1),n*N,n*N)*J0 ...
    - sparse(1:n*N,1:n*N,repmat(Q(1:N),n,1),n*N,n*N) ...
    - U*DQ(1:N,:) - J0*U*DQ(1+(n-1)*N:n*N,:);

J = Lap + p(2)*DQ + p(3)*DC;
end

end