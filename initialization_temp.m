function [init] =initialization_temp(N,sqn,x)

temp=zeros(N,1);
    for i=1:(sqn*sqn) 
        temp(i)=5*sin((x(1,i)+x(2,i))*100);%fix((i-1)/sqn)*v_0;
    end
 
init=temp;