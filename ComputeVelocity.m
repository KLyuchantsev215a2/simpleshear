function [velocity]=ComputeVelocity(dt,v,SIG,nabla_W,V,N,m,viscosity,rho_0)
velocity=v;
nabla_Wij=[0,0];
for i=1:N
    for j=1:N
        for beta=1:2
            for alpha=1:2
                
                nabla_Wij=nabla_W(1:2,i,j);
               % viscosity=dt*eps*h*cs*V(j)*(v(alpha,i))*(Hessian_W_cor(1,i,j)+Hessian_W_cor(2,i,j));   % - v(alpha,j)
velocity(alpha,i)=velocity(alpha,i)+dt*(((m/(rho_0*rho_0))*SIG(alpha,beta,j))*nabla_Wij(beta));
             end
        end
    end
    
    for alpha=1:2
         velocity(alpha,i) = velocity(alpha,i) + dt*viscosity(alpha,i);
    end
   
 end
