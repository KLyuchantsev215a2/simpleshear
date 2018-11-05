function [Fn1]=ComputeF(L,N,dt,F)

 for i = 1:N    
         dtLL = dt* L(1:2,1:2,i);
         F(1:2,1:2,i)= expm(dtLL)*F(1:2,1:2,i);
end

Fn1=F;