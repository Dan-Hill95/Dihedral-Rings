%Dan J Hill (2021) - Simple D_{2k} patch matching condition:
%N is the truncation order/dimension index
%k is the lattice index
%a is an (N+1) vector; we will solve for this close to some initial guess
function [F,J]=match(a,k)
%Defining N and F
N=length(a)-1;
F=sparse(1,N+1);
for i=1:N+1
   F(i)= a(i);
   for j=2:N+2-i
   F(i)= F(i)-2*cos((k*pi/3)*(i-j))*a(j)*a(i+j-1);    
   end
   for j=1:i
   F(i)= F(i)-cos((k*pi/3)*(i-2*j+1))*a(j)*a(i-j+1);        
   end
end
end