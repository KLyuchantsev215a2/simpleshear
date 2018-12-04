function [acc,F,SIG]=ComputeAcceleration(x,v,V,N,m,eps1,eps2,h,cs_0,nabla_W_cor_0,X_old,mu,k,dh)

F=ComputeF(V,x,nabla_W_cor_0,N,X_old);  
SIG=ComputeStress(F,mu,k,N);
acc=zeros(2,N);

[W_cor,nabla_W_cor,Hessian_W_cor]=ComputeW_final(x,V,N,h,dh);
V=computeV(N,W_cor,m);  

nabla_Wij=[0,0];

viscosity=ComputeViscocity(v,V,eps1,h,Hessian_W_cor,cs_0,N,eps2);

for i=1:N
    for j=1:N
        for alpha=1:2 
            for beta=1:2                   
                nabla_Wij=nabla_W_cor(1:2,i,j);
                acc(alpha,i)=acc(alpha,i)+((V(j)))*SIG(alpha,beta,j)*nabla_Wij(beta);
             end
        end
    end
   
    for alpha=1:2
         acc(alpha,i)=acc(alpha,i)*V(i)/m;
         acc(alpha,i) = acc(alpha,i) + viscosity(alpha,i);
    end
end
