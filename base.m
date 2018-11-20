clear;    

rho_0 =120;
v_0 = 0.1;
Time = 10;
sqn=5;
l=0.01;
N=sqn*sqn;
S=l*l;
m=rho_0*S/N;

nu=0.4;
mu = 10;  
k=2*mu*(1+nu)/(3*(1-2*nu));
E=9*k*mu/(3*k+mu);   % модуль Юнга

cs_0=sqrt((E+4/3*mu)/rho_0);

h=1*(m/rho_0)^(1/2);%k увеличен
dt=0.0005;
dh=0.0000001;
eps1=1/4;%-100;
eps2=1/8;%-50;%1/5;

V=m/rho_0*ones(N,1);%m/rho_0;
x=initialization_x(N,sqn,l);    
v=initialization_v(N,sqn,v_0,x);
viscosity=zeros(2,N);
X_old=x;

F=zeros(2,2,N);
SIG=zeros(2,2,N);


for i = 1:N
    F(1:2,1:2,i)=eye(2);
end

SIG=ComputeStress(F,mu,k,N);
W_cor=zeros(N,N);
nabla_W_cor=zeros(2,N,N);
Hessian_W_cor=zeros(2,N,N);

[W_cor,nabla_W_cor_0,Hessian_W_cor]=ComputeW_final(x,V,N,h,dh);
for n = 1:fix(Time/dt)
    x_old=x;
    v_old=v;
    [acc_old,F_old,SIG_old]=ComputeAcceleration(x_old,v_old,V,N,m,eps1,eps2,h,cs_0,nabla_W_cor_0,X_old,mu,k,dh);
    x_star = x_old + dt*v_old;
    v_star = v_old + dt*acc_old;
    [acc_star,F,SIG] = ComputeAcceleration(x_star,v_star,V,N,m,eps1,eps2,h,cs_0,nabla_W_cor_0,X_old,mu,k,dh);
    x_star_star = x_star + dt*v_star;
    v_star_star = v_star + dt*acc_star;
    x = 0.5*(x_old + x_star_star);
    v = 0.5*(v_old + v_star_star);
    plotmy=myplot(x,V,F,N,SIG,l,v);
    life_time=n*dt;
end
