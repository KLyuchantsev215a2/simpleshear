 function [V] = computeV(N,W_cor,m)

V=zeros(N,1);
rho=zeros(N,1);
for i = 1:N
    for j = 1:N
       rho(i)=rho(i)+m*W_cor(i,j); 
    end
end

for i = 1:N
<<<<<<< HEAD
    V(i,1)=m/V(i,1);
=======
    V(i)=m/rho(i);
>>>>>>> 33acd43b5d49aeb99e8f8d2532a5815042314075
end
