%% Direct Frequency Domain Analysis 
% Running on the unit square with homogenous discretization
% -cΔu + c*u = f 
% (L - w*w*M) + i(mu*L + 10*w*I)) u = f 
% Both sides are multiplied with the step-length

function [C, c] = DFDA(m,omega,mu) 
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

%mu = 0.02;
%omega = pi; 
disp('note, omega and mu are fixed')

C = h*h*(K - omega*omega*mu*Imm + 1i*((omega*10*Imm)+mu*K));

% J = [1:(m*m)]';
% c = h*h*(1-1i)*J./(m*m)./(h*(J+1).^2); 
% c = h*h*(1+1i)*(C*ones(size(C(:,1)))); 
c =  h*h*((1-2*rand(m*m,1))+1i*(1-2*rand(m*m,1))); 
