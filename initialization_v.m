function [init] =initialization_v(N,sqn,v_0,x)

v=zeros(2,N);
    for i=1:(sqn*sqn) 
        v(1,i)=0;%-x(2,i)*v_0;
        v(2,i)=0;
    end
    

init=v;