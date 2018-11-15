clear;    

rho_0 =0.012;
v_0 = 0.1;
Time = 10;
sqn=8;
l=0.01;
N=sqn*sqn;
S=l*l;
m=rho_0*S/N;

nu=0.4;
mu = 10;  
k=2*mu*(1+nu)/(3*(1-2*nu));
E=9*k*mu/(3*k+mu);   % ������ ����

cs_0=sqrt((E+4/3*mu)/rho_0);

h=4*(m/rho_0)^(1/2);%k ��������
dt=0.000001;
dh=0.0000001;
eps=0.0025;


V=m/rho_0*ones(N,1);%m/rho_0;
x=initialization_x(N,sqn,l);    
v=initialization_v(N,sqn,v_0,x);
temp=initialization_temp(N,sqn,x);
viscosity=zeros(2,N);
X_old=x;

F=zeros(2,2,N);
SIG=zeros(2,2,N);
P=zeros(2,2,N);%�����_�������� 1
f=zeros(2,N);%������ ���������� ���

for i = 1:N
    F(1:2,1:2,i)=eye(2);
end

SIG=ComputeStress(F,mu,k,N);

    
for n = 1:fix(Time/dt)
    
    W_cor=zeros(N,N);
    W_cor_1per=zeros(N,N);
    W_cor_2per=zeros(N,N);
    nabla_W_cor=zeros(2,N,N);
    Hessian_W_cor=zeros(2,N,N);
    W_cor=ComputeW_cor(N,x,x,V,h);
    
    x_per1=x;
    x_per2=x;
    x_per1_inv=x;
    x_per2_inv=x;
    
    x_per1(1,1:N)=x(1,1:N)+dh;
    x_per2(2,1:N)=x(2,1:N)+dh;
    x_per1_inv(1,1:N)=x(1,1:N)-dh;
    x_per2_inv(2,1:N)=x(2,1:N)-dh;
    
    W_cor_1per=ComputeW_cor(N,x,x_per1,V,h);
    W_cor_2per=ComputeW_cor(N,x,x_per2,V,h);
    W_cor_1per_inv=ComputeW_cor(N,x,x_per1_inv,V,h);
    W_cor_2per_inv=ComputeW_cor(N,x,x_per2_inv,V,h);
    
    nabla_W_cor(1,1:N,1:N)=(W_cor_1per-W_cor)/dh;
    nabla_W_cor(2,1:N,1:N)=(W_cor_2per-W_cor)/dh;
    
    Hessian_W_cor(1,1:N,1:N)=(W_cor_1per-2*W_cor+W_cor_1per_inv)/(dh*dh);
    Hessian_W_cor(2,1:N,1:N)=(W_cor_2per-2*W_cor+W_cor_2per_inv)/(dh*dh);
    
    if (n==1)
        nabla_W_cor_0=nabla_W_cor;
    end
    
    L=ComputeL(v,V,nabla_W_cor,N);
   % temp=ComputeTemp(temp,V,eps,Hessian_W_cor,N);
    V=computeV(N,W_cor,m); 
    viscosity=ComputeViscocity(v,V,eps,h,Hessian_W_cor,cs_0,N);
    v=ComputeVelocity(dt,v,SIG,nabla_W_cor,V,N,m,viscosity);

    for i = 1:N
        x(1,i)=x(1,i)+dt*v(1,i);
        x(2,i)=x(2,i)+dt*v(2,i);
    end
    
    F=ComputeF(V,x,nabla_W_cor_0,N,X_old);  
    SIG=ComputeStress(F,mu,k,N);
    P=ComputeKirchhoff(F,SIG,N);
    f=ComputeForse(V,P,nabla_W_cor,N);
    plotmy=myplot(x,V,F,N,SIG,l,temp);
%    fsum=[0;0];
%     for i=1:N
%         fsum(1:2)=f(1:2,i)+fsum(1:2);
%     end
%     fsum=fsum;
end
