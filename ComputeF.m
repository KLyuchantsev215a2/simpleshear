function [Fn1]=ComputeF(V,x,nabla_W0,N)
Fn1=zeros(2,2,N);
for i=1:N
    for j=1:N
        for beta=1:2
            for alpha=1:2
            Fn1(alpha,beta,i)= Fn1(alpha,beta,i)+(V(j)* (x(alpha,j)-x(alpha,i)) *nabla_W0(beta,i,j));  
            end
        end
    end
end
    
             

