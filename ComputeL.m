function [Ln1]=ComputeL(v,V,nabla_W_cor,N)

L=zeros(2,2,N);
nabla_W_tmp=[0,0];
L_tmp=zeros(2,2);
v_tmp=[0,0];


for i = 1:N
    v_tmpi=v(1:2,i);
    for j = 1:N
         
         v_tmp=v(1:2,j);
         
         nabla_W_tmp=nabla_W_cor(1:2,i,j); 
        
         for beta=1:2            
            L_tmp(1,beta)=L_tmp(1,beta)+V(1,j)*(v_tmp(1)-v_tmpi(1))*nabla_W_tmp(beta);  
            L_tmp(2,beta)=L_tmp(2,beta)+V(1,j)*(v_tmp(2)-v_tmpi(2))*nabla_W_tmp(beta); 
         end
        
    end
    L(1:2,1:2,i)=L_tmp;
    L_tmp=zeros(2,2);
end

Ln1=L;