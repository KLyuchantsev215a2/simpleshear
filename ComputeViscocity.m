function [viscosity]=ComputeViscocity(v,V,eps,h,Hessian_W_cor,cs,N)
viscosity=zeros(2,N);
for i=1:N
    for j=1:N
        for alpha=1:2
            viscosity(alpha,i)=viscosity(alpha,i)+eps*h*cs*V(j)*(v(alpha,j)-v(alpha,i))*(Hessian_W_cor(1,i,j)+Hessian_W_cor(2,i,j));
       end
    end
 end

