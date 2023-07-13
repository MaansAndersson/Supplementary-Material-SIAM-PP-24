%% Padé Bai
% Running on the unit square with homogenous discretization
% -cΔu + c*u = f 
% Both sides are multiplied with the step-length

function [C, c] = PadeBai(m,sigma0,sigma1) 
% Dv = -[-1/12*e 4/3*e -5/2*e 4/3*e -1/12*e];
%L  = spdiags(Dv,-2:2,m,m)/(h*h);

h = 1/(m-1);

e  = ones(m,1);
Dv = -[e -2*e e];
L  = spdiags(Dv,-1:1,m,m)/(h*h);

Imm = speye(m*m);
Im = speye(m);

% on a 2D grid
K = kron(Im,L) + kron(L,Im);

C = h*h*(K + (3-sqrt(3))/h*Imm + 1i*(K +(3+sqrt(3))/h*Imm)); 

J = [1:(m*m)]';
c = h*h*(1-1i)*J./(m*m)./(h*(J+1).^2); 
% c = h*h*(1+1i)*(C*ones(size(C(:,1)))); 



