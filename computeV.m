 function [V] = computeV(N,W_cor,m)

V=zeros(N,1);
for i = 1:N
    for j = 1:N
       V(i)=V(i)+W_cor(i,j); 
    end
end

for i = 1:N
    V(i,1)=m/V(i,1);
end
