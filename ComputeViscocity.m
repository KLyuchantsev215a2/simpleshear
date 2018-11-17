function [viscosity]=ComputeViscocity(v,V,eps1,h,Hessian_W_cor,cs,N,eps2)
viscosity_1=zeros(2,N);
viscosity_2=zeros(2,N);
viscosity=zeros(2,N);
viscosity_22=zeros(2,N);
for i=1:N
    for j=1:N
        for alpha=1:2
            viscosity_22(alpha,i)=viscosity_22(alpha,i)+eps1*h*cs*V(j)*(v(alpha,j)-v(alpha,i))*(Hessian_W_cor(1,i,j)+Hessian_W_cor(2,i,j));
       end
    end
end
 
 for j=1:N
    for i=1:N
        for alpha=1:2
            viscosity_1(alpha,j)=viscosity_1(alpha,j)+V(i)*(v(alpha,i)-v(alpha,j))*(Hessian_W_cor(1,j,i)+Hessian_W_cor(2,j,i));
       end
    end
 end

for i=1:N
    for j=1:N
        for alpha=1:2
            viscosity_2(alpha,i)=viscosity_2(alpha,i)+V(j)*(v(alpha,j)-v(alpha,i))*(Hessian_W_cor(1,i,j)+Hessian_W_cor(2,i,j));
       end
    end
end

for i=1:N
    for j=1:N
        for alpha=1:2
            viscosity(alpha,i)=-eps2*h*h*h*cs*V(j)*(viscosity_1(alpha,j)- viscosity_2(alpha,i))*(Hessian_W_cor(1,i,j)+Hessian_W_cor(2,i,j))+viscosity_22(alpha,i);
       end
    end
end


