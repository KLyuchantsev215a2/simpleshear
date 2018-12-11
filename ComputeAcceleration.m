function [acc,F,SIG]=ComputeAcceleration(x,v,V,N,m,eps1,eps2,h,cs_0,nabla_W_cor_0,X_old,mu,k,dh,rho_0)

F=ComputeF(V,x,nabla_W_cor_0,N,X_old);  
SIG=ComputeStress(F,mu,k,N);
acc=zeros(2,N);

[W_cor,nabla_W_cor,Hessian_W_cor]=ComputeW_final(x,V,N,h,dh);
V=computeV(N,W_cor,m);  

SIG_tmp=zeros(2);
viscosity=ComputeViscocity(v,V,eps1,h,Hessian_W_cor,cs_0,N,eps2);

for i=1:N
    for j=1:N
            for beta=1:2    
                 for alpha=1:2 
                acc(alpha,i)=acc(alpha,i)-((V(j)))*SIG(alpha,beta,j)*nabla_W_cor(beta,j,i);
             end
        end
    end
   
    for alpha=1:2
         acc(alpha,i)=acc(alpha,i)/rho_0;
         acc(alpha,i) = acc(alpha,i) + viscosity(alpha,i);
    end
end
