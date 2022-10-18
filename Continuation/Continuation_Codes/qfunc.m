% solve A_{ss} = -A_s/s + A/4s^2 + c0*A + c3*A^3
% solved using 2nd order finite differences
function A=qfunc(c0, c3, Li,Ls,Ns)
%c0 = 1;     % c0
%c3 = -1;    % c3

% Li = 0;
% Ns = 1000;   % number of s discretisation points
%Ls = 10;    % compute on s in [0,Ls]
s  = linspace(Li,Ls,Ns)'; % s mesh
si = s(2:Ns-1); % interior s points
ds = s(2)-s(1); % stepsize of mesh

A0 = 2*sqrt(s(2:Ns-1)).*sech(s(2:Ns-1)); % initial conditon for Newton

e = ones(Ns-2,1);                       % differentiation matrices on interior mesh si
Ds  = spdiags([-e 0*e e], -1:1, Ns-2, Ns-2)/(2*ds);
Dss = spdiags([e -2*e e], -1:1, Ns-2, Ns-2)/ds^2;

% figure;
% plot(si,A0,'b',si,Ds*A0,'r',si,Dss*A0,'m');
% plot(si,Ds*A,'b',si,(sech(si)./(2*sqrt(si)) - 2*sqrt(si).*sech(si).*tanh(si) ),'r');
% plot(si,Dss*A,'b',si,((4*si.^2-1).*cosh(si).^2-4*sinh(si).*si.*cosh(si)-8*si.^2)./(2*si.^(3/2).*cosh(si).^3),'r')
% plot(si,Ds*sin(si*pi/Ls),'b',si,pi/Ls*cos(si*pi/Ls),'--r');

options = optimoptions('fsolve','Display','iter');
[Aj,fval,exitflag,output,jacobian] = fsolve(@(A) radial_GZ_eqn(A,c0,c3,Ds,Dss,si),A0,options);
A = [Aj(1);Aj;0];
end