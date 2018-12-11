function [U] = ComputeEnergy(Fp,mu,N)
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
    J = det(F);
    B=F*F';
    Ciso = J^(-2/3)*B;    % isochoric part of  right Cauchy–Green deformation tensor
    U(i)=mu/2*(trace(Ciso)-3);   
end