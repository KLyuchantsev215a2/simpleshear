clear;    
<<<<<<< HEAD
rho_0 =1.2;
v_0 = 10;
=======
rho_0 =0.12;
v_0 = 0.1;
>>>>>>> 33acd43b5d49aeb99e8f8d2532a5815042314075
Time = 10;
sqn=6;
l=0.01;
N=sqn*sqn;
S=l*l;

nu=0.35;
mu=10;
k=2*mu*(1+nu)/(3*(1-2*nu));


m=rho_0*S/N;
<<<<<<< HEAD
h=4.8*(m/rho_0)^(1/2);
dt=0.00001;
=======
h=1*(m/rho_0)^(1/2);
dt=0.000005;
>>>>>>> 33acd43b5d49aeb99e8f8d2532a5815042314075
dh=0.00000001;
V=m/rho_0*ones(N,1);%m/rho_0;
x=initialization_x(N,sqn,l);    
v=initialization_v(N,sqn,v_0,x);

F=zeros(2,2,N);
SIG=zeros(2,2,N);
           
for i = 1:N
    F(1:2,1:2,i)=eye(2);
    SIG(1:2,1:2,i)=ComputeStress(F(1:2,1:2,i),mu,k);
end

    
for n = 1:fix(Time/dt)
    W_cor=zeros(N,N);
    W_cor_1per=zeros(N,N);
    W_cor_2per=zeros(N,N);
    nabla_W_cor=zeros(2,N,N);
    W_cor=ComputeW_cor(N,x,x,V,h);
    x_per1=x;
    x_per2=x;
    
    x_per1(1,1:N)=x_per1(1,1:N)+dh;
    x_per2(2,1:N)=x_per2(2,1:N)+dh;
    
    W_cor_1per=ComputeW_cor(N,x,x_per1,V,h);
    W_cor_2per=ComputeW_cor(N,x,x_per2,V,h);
    
    nabla_W_cor(1,1:N,1:N)=(W_cor_1per-W_cor)/dh;
    nabla_W_cor(2,1:N,1:N)=(W_cor_2per-W_cor)/dh;
    
    if (n==1)
        nabla_W_cor_0=nabla_W_cor;
    end
    
    L=ComputeL(v,V,nabla_W_cor,N);
        
<<<<<<< HEAD
 
    V=computeV(N,W_cor,m); 
    
       for i = 1:N
           v(2,i)=(x(2,i))*v_0;
       end
        
        for i = 1:N
            x(2,i)=x(2,i)+dt*v(2,i);
        end
    
        
=======
   % v=ComputeVelocity(dt,v,SIG,nabla_W_cor,V,m,N);
    
    for i = 1:N
        x(1,i)=x(1,i)+dt*v(1,i);
        x(2,i)=x(2,i)+dt*v(2,i);
    end
    
     V=computeV(N,W_cor,m);
>>>>>>> 33acd43b5d49aeb99e8f8d2532a5815042314075
    
%      for i = 1:N   
%         dtLL = dt* L(1:2,1:2,i);
%         F(1:2,1:2,i)= expm(dtLL)*F(1:2,1:2,i);
%      end
    F=ComputeF(V,x,nabla_W_cor_0,N);
     
    for i = 1:N   
          SIG(1:2,1:2,i)=ComputeStress(F(1:2,1:2,i),mu,k);
    end
      
<<<<<<< HEAD
    plotmy=myplot(x,V,F,N);
   
=======
    plotmy=myplot(x,V,F,N,SIG);
>>>>>>> 33acd43b5d49aeb99e8f8d2532a5815042314075
    

end
