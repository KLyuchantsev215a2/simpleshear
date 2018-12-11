function [temperature]=ComputeTemp(temp,V,eps,Hessian_W_cor,N)
temperature = temp;
for i=1:N
    for j=1:N        
            temperature(i)=temperature(i)+eps*V(j)*(temp(j)-temp(i))*(Hessian_W_cor(1,i,j)+Hessian_W_cor(2,i,j));       
    end
 end

