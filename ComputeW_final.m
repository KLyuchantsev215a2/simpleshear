function [W_cor,nabla_W_cor, Hessian_W_cor]=ComputeW_final(x,V,N,h,dh)

        W_cor=zeros(N,N);
        W_cor_1per=zeros(N,N);
        W_cor_2per=zeros(N,N);
        nabla_W_cor=zeros(2,N,N);
        Hessian_W_cor=zeros(2,N,N);
    
    
        x_per1=x;
        x_per2=x;
        x_per1_inv=x;
        x_per2_inv=x;
    
        x_per1(1,1:N)=x_per1(1,1:N)+dh;
        x_per2(2,1:N)=x_per2(2,1:N)+dh;
        x_per1_inv(1,1:N)=x_per1_inv(1,1:N)-dh;
        x_per2_inv(2,1:N)=x_per2_inv(2,1:N)-dh;
    
        W_cor=ComputeW_cor(N,x,x,V,h);
    
        W_cor_1per=ComputeW_cor(N,x,x_per1,V,h);
        W_cor_2per=ComputeW_cor(N,x,x_per2,V,h);
        W_cor_1per_inv=ComputeW_cor(N,x,x_per1_inv,V,h);
        W_cor_2per_inv=ComputeW_cor(N,x,x_per2_inv,V,h);
    
        nabla_W_cor(1,1:N,1:N)=(W_cor_1per-W_cor)/dh;
        nabla_W_cor(2,1:N,1:N)=(W_cor_2per-W_cor)/dh;
        
        Hessian_W_cor(1,1:N,1:N)=(W_cor_1per-2*W_cor+W_cor_1per_inv)/(dh*dh);
        Hessian_W_cor(2,1:N,1:N)=(W_cor_2per-2*W_cor+W_cor_2per_inv)/(dh*dh);

