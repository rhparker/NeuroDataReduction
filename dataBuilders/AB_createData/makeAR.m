function [S] = makeAR(nC,T,lam)
% Create inhomogeneous Poisson processes using an autoregressive 
% process as described in RD 2017
%
% INPUTS
%    nC     = # cells
%    T      = length of time
%    lam    = mean rate
%nC  = 10;

D       = 0.9*eye(nC);      % Coefficient matrix
sigS    = 0.01;             % Std dev of noise 
%T       = 1400;             % Final time

% For sigmoidal nonlinearity
%lam     = 5;                % Mean rate
v       = 0.2;              % Slope for sigmoid
sbar    = 0.02;              % Subtract

S       = zeros(nC,T+1);
S(:,1)  = sigS*randn(nC,1);   %Initial condition

for j=1:T
    S(:,j+1)=D*S(:,j) + sigS*randn(nC,1);
end

S = lam*(1+erf(v*(S-sbar)/sigS));

end

