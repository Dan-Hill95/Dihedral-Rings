% Spatial coordinates: r direction
N = 1000; % Number of mesh points, must be even! 
L = 100; % domain truncation
h = L/(N-1); r = (0:N-1)'*h;

%% Setup differentiation matrices

% compute radial Swift-Hohenberg operator differentiation matrix
e = ones(N,1);

% d_r
D = sparse(1:N-1,[2:N-1 N],ones(N-1,1)/2,N,N); 
D = (D - D')/(h);
D(1,2) = 0;
D(N,N-1) = 0;

% d_rr
D2 = sparse(1:N-1,[2:N-1 N],ones(N-1,1),N,N) - sparse(1:N,[1:N],e,N,N);
D2 = (D2 + D2');

D2(1,2)=2;D2(N,N-1)=2;

D2 = D2/h^2;

% 1/r
r(1) = 1;
R = sparse(1:N,[1:N],1./r,N,N);
r(1) = 0;
e = sparse(1:N,[1:N],e,N,N); % I

LN = D2 + R*D + e; % d_rr + d_r/r + I
LN(1,1) = 2*D2(1,1) + 1; LN(1,2) = 2*D2(1,2); % 2d_rr at r= 0
LN = -LN^2; % -(1+laplacian)^2

LM=sparse((n-1)*N,(n-1)*N);

for i=1:n-1
LM(1+(i-1)*N:i*N,1+(i-1)*N:i*N) = D2 + R*D + e - (i*m)^2*(R.^2).*e; % d_rr + d_r/r + I-(i*m/r)^2
LM(1+(i-1)*N,1+(i-1)*N) = (2-((m*i)^2)/2)*D2(1,1) + 1; LM(1+(i-1)*N,2+(i-1)*N) = (2-((m*i)^2)/2)*D2(1,2); % 2d_rr at r= 0
LM(1+(i-1)*N:i*N,1+(i-1)*N:i*N) = -LM(1+(i-1)*N:i*N,1+(i-1)*N:i*N)^2; % -(1+laplacian)^2
end

mesh_params.N  = N;

mesh_params.L  = L; mesh_params.m = m; mesh_params.e = e;

mesh_params.r  = r; mesh_params.R = R;

mesh_params.D  = D; mesh_params.D2 = D2; mesh_params.LN = LN; mesh_params.LM = LM; 
