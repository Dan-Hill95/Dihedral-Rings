function F = radial_GZ_eqn(A,c0,c3,Ds,Dss,s)

F = Dss*A + (Ds*A)./s - A./(4*s.^2) - c0*A - c3*A.^3;
