 function [V] = computeV(N,W_cor)

V=zeros(N,1);
for i = 1:N
    for j = 1:N
       V(i)=V(i)+W_cor(i,j); 
    end
end

