 function [V] = computeV(N,W_cor)

V=zeros(1,N);
for i = 1:N
    for j = 1:N
       V(1,i)=V(1,i)+W_cor(i,j); 
    end
end

