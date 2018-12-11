clear;    

% fig = figure();
% % создание первого пустого кадра
% set(fig,'Position',[350,200,700,700]);
% frame = getframe(fig);
% [im,map] = rgb2ind(frame.cdata,4);
% imwrite(im,map,'animation3.gif','DelayTime',0,'Loopcount',inf);

rho_0 =3;
v_0 = 1;
Time = 10;
sqn=5;
l=1;
N=sqn*sqn;
S=l*l;
m=rho_0*S/N;

nu=0.4;
mu = 10;  
k=2*mu*(1+nu)/(3*(1-2*nu));
E=9*k*mu/(3*k+mu);   % модуль Юнга

cs_0=sqrt((E+4/3*mu)/rho_0);

h=1.4*(m/rho_0)^(1/2);%k увеличен
dt=0.1*h/(cs_0+v_0);
dh=0.0000001;
eps1=0;%-100;
eps2=0;%-50;%1/5;

V=m/rho_0*ones(N,1);%m/rho_0;
V_0=m/rho_0*ones(N,1);%m/rho_0;
x=initialization_x(N,sqn,l);    
v=initialization_v(N,sqn,v_0,x,l);
viscosity=zeros(2,N);
X_old=x;

F=zeros(2,2,N);
L=zeros(2,2,N);
SIG=zeros(2,2,N);
U=zeros(N,1);
Energy_time=zeros(fix(Time/dt),1);
U_time=zeros(fix(Time/dt),1);%potential energy
T_time=zeros(fix(Time/dt),1);%kinetic energy
time=zeros(fix(Time/dt),1);

for i = 1:N
    F(1:2,1:2,i)=eye(2);
end

SIG=ComputeStress(F,mu,k,N);
W_cor=zeros(N,N);
nabla_W_cor=zeros(2,N,N);
Hessian_W_cor=zeros(2,N,N);

[W_cor,nabla_W_cor_0,Hessian_W_cor]=ComputeW_final(x,V,N,h,dh);
for n = 1:fix(Time/dt)
    L=zeros(2,2,N);
%     if(fix(n/200)==n/200)
%         X_old=x;
%         [W_cor,nabla_W_cor_0,Hessian_W_cor]=ComputeW_final(x,V,N,h,dh);
%     end
%     x_old=x;
%     v_old=v;
%     [acc_old,F_old,SIG_old]=ComputeAcceleration(x_old,v_old,V,N,m,eps1,eps2,h,cs_0,nabla_W_cor_0,X_old,mu,k,dh);
%     x_star = x_old + dt*v_old;
%     v_star = v_old + dt*acc_old;
%     [acc_star,F,SIG] = ComputeAcceleration(x_star,v_star,V,N,m,eps1,eps2,h,cs_0,nabla_W_cor_0,X_old,mu,k,dh);
%     x_star_star = x_star + dt*v_star;
%     v_star_star = v_star + dt*acc_star;
%     x = 0.5*(x_old + x_star_star);
%     v = 0.5*(v_old + v_star_star);
%    [W_cor,nabla_W_cor,Hessian_W_cor]=ComputeW_final(x,V,N,h,dh);
%    L=ComputeL(v,V,nabla_W_cor,N);
      x_old=x;
      v_old=v;
      [acc_old,F_old,SIG_old]=ComputeAcceleration(x_old,v_old,V,N,m,eps1,eps2,h,cs_0,nabla_W_cor_0,X_old,mu,k,dh,rho_0);
      x_n_1=x_old+dt*v_old;
      v_n_1=v_old+dt*acc_old;
      [acc_n_1,F,SIG]=ComputeAcceleration(x_n_1,v_n_1,V,N,m,eps1,eps2,h,cs_0,nabla_W_cor_0,X_old,mu,k,dh,rho_0);
      x_n_2=x_n_1+dt*v_n_1;
      v_n_2=v_n_1+dt*acc_n_1;
      x_n_1_2=3/4*x_old+1/4*x_n_2;
      v_n_1_2=3/4*v_old+1/4*v_n_2;
      [acc_n_1_2,F,SIG]=ComputeAcceleration(x_n_1_2,v_n_1_2,V,N,m,eps1,eps2,h,cs_0,nabla_W_cor_0,X_old,mu,k,dh,rho_0);
      x_n_3_2=x_n_1_2+dt*v_n_1_2;
      v_n_3_2=v_n_1_2+dt*acc_n_1_2;
      x=1/3*x_old+2/3*x_n_3_2;
      v=1/3*v_old+2/3*v_n_3_2;
      U=ComputeEnergy(F,mu,k,N);
      Energy=0;
      U_energy = 0;
      T_energy = 0;
      for i=1:N
          % Energy=Energy+U(i)*V_0(i)+1/2*rho_0*V_0(i)*dot(v(1:2,i),v(1:2,i));%+U(i)*V(i)
          U_energy = U_energy + U(i)*V_0(i);
          %T_energy = T_energy + 1/2*rho_0*V_0(i)*dot(v(1:2,i),v(1:2,i)); 
          T_energy = T_energy + 1/2*rho_0*V_0(i)*(v(1,i)^2+ v(2,i)^2);
      end
      U_time(n)=U_energy;%potential energy
      T_time(n)= T_energy;
      Energy_time(n)=U_energy + T_energy;
      time(n)=n*dt;
%         subplot(2,2,1);
%         x_coord =time;
%         y_coord = U_time;
%         plot(x_coord,y_coord);
%         subplot(2,2,2);
%         x_coord =time;
%         y_coord = T_time;
%         plot(x_coord,y_coord);
%           subplot(2,2,3);
%         x_coord =time;
%         y_coord = Energy_time;
%         plot(x_coord,y_coord);
%          pause(0.0000001);
     plotmy=myplot(x,V,F,N,SIG,l,v,Energy_time,time);%%n,im,frame,map,fig);
      life_time=n*dt;
end
% x_coord =time;
% y_coord = Energy_time;
% subplot(2,2,1);
% plot(x_coord,y_coord);

