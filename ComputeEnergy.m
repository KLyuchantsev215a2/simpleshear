function [U] = ComputeEnergy(Fp,mu,k,N)
% Subroutine which computes Cauchy stress tensor as a function
% of the deformation gradient F
% input:  F = deformation gradient
%         mu = shear modulus
%          k = bulk modulus
% output: stress = Cauchy stresss (true stress) 
%
%   neo-Hookean material
%
U=zeros(N,1);
for i = 1:N        
    F = Fp(1:2,1:2,i);
    F(3,3) = 1;
    
   % F_0 = [1,0.1,0;0,1,0;0,0,1];  % new unloaded state
   % F = F*F_0^(-1);               % new unloaded state
    
    J = det(F);
    B=F*F';
    Biso = J^(-2/3)*B;    % isochoric part of  right Cauchy–Green deformation tensor
    U(i)=mu/2*(trace(Biso)-3)+k/50*(J^5+J^(-5)-2);   
end