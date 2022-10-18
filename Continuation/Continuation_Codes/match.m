%Dan J Hill (2022) - D_{k} ring matching condition:
%N is the truncation order/dimension index
%m is the lattice index
%a is an (N+1) vector; we will solve for this close to some initial guess
function [F,J]=match(a,m)
%Defining N and F
N=length(a)-1;
F=a;
% if N==1
%     F(1)=a(1) - a(1)^3 - 2*(2+(-1)^m)*a(1)*a(2)^2;
%     F(2)=a(2) - 3*a(2)^3 - (2+(-1)^m)*a(1)^2*a(2);
% else
for n=0:N
   for i=-N:N
       for j= -N:N-abs(i-n)
           %F_n = a_n - sum_i=-N^N sum_j=-N^N-|n-i| (-1)^[(m/2)*(|i| + |j| - |j+|n-i||)]*a_{|i|}*a_{|j|}*a_{|j+|n-i||}
       F(n+1) = F(n+1) - ((-1)^((m/2)*(abs(i)+abs(j) - abs(j + abs(n-i))-n)))*a(abs(i)+1)*a(abs(j)+1)*a(abs(j+abs(n-i))+1);
% for i = 0:n
% for j = 1:N-i
% q1()=
% end
% for j = 1:N-i
% q2()=
% end
% F(n+1)=F(n+1)- \sum_{i=0}^{n} a_{n-i} \left\{\sum_{j= 1}^{N-i} \big[1 + (-1)^{-m i}\big]  a_{j}a_{j+i} + \sum_{j= 0}^{i} (-1)^{-m j} a_{j}a_{i-j}\right\}\\ &\quad + \sum_{i=1}^{N-n}  a_{i+n}\left\{\sum_{j= 1}^{N-i} \big[1+(-1)^{mi} \big] a_{j}a_{j+i}+ \sum_{j= 0}^{i} (-1)^{m j} a_{j}a_{i-j}\right\}\\ &\quad + (-1)^{m n}\sum_{i=n+1}^{N} a_{i-n} \left\{\sum_{j= 1}^{N-i} \big[1 + (-1)^{m i}\big] a_{j}a_{j+i} + \sum_{j= 0}^{i} (-1)^{m j} a_{j}a_{i-j}\right\}\\ &\quad + (-1)^{m n} \sum_{i=1}^{n} a_{i-n+N} \sum_{j= i}^{N} (-1)^{m j} a_{j}a_{i+N-j},\\;


       end
   end
end
if nargout > 1 %Calculating the Jacobian
   J = speye(N+1);
   delta = speye(N+1);
   %sparse(1+j*N:n*N,1:(n-j)*N,repmat(Q(1+j*N:(j+1)*N),n-j,1),n*N,n*N)
for n=0:N
   for i=-N:N
       for j= -N:N-abs(i-n)
           for k=1:N+1
           %F_n = a_n - sum_i=-N^N sum_j=-N^N-|n-i| (-1)^[(m/2)*(|i| + |j| - |j+|n-i||)]*a_{|i|}*a_{|j|}*a_{|j+|n-i||}
       J(n+1,k) = J(n+1,k) - ((-1)^((m/2)*(abs(i)+abs(j) - abs(j + abs(n-i))-n)))*delta(abs(i)+1,k)*a(abs(j)+1)*a(abs(j+abs(n-i))+1) - ((-1)^((m/2)*(abs(i)+abs(j) - abs(j + abs(n-i))-n)))*a(abs(i)+1)*delta(abs(j)+1,k)*a(abs(j+abs(n-i))+1) - ((-1)^((m/2)*(abs(i)+abs(j) - abs(j + abs(n-i))-n)))*a(abs(i)+1)*a(abs(j)+1)*delta(abs(j+abs(n-i))+1,k);
           end
       end
   end
end
end
end
