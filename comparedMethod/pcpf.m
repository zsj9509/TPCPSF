function [L, S] = pcpf(M, W, X, Y, varargin)
%Runs the principal component pursuit using side information with features
%algorithm implemented using inexact augmented Lagrange multiplier method

%%Initialisation
[m, n] = size(M);
d = [size(X,2),size(Y,2)];
M_fnorm = norm(M,'fro');

if nargin > 5
    lambda = varargin{1};
    kappa = varargin{2};
else
    lambda = 1/sqrt(max(m,n));
    
    if nargin > 4
        kappa = varargin{1};
    else
        kappa = 0;
    end
end

mu = 1/lansvd(M, 1, 'L');
alpha = 1.1;
mu_max = 1e18;

Z = zeros(m, n);

H = zeros(d);
E = H;
N = H;

tol = 1e-7; % pcpsf的tol
% tol = 1e-8;   % tprca的tol

iter = 0;
iter_max = 1000;

hasConverged = false;

%Calculation
while ~hasConverged && iter < iter_max
    iter = iter + 1;
    
    S = get_S_update;
    H = get_H_update;
    E = get_E_update;
    
    %cache variables
    Temp1 = M - S - X*H*Y';
    Temp2 = H - E - X'*W*Y;
    
    Z = get_Z_update(Temp1);
    N = get_N_update(Temp2);
    
    mu = min(mu*alpha, mu_max);  
    % 这个是PSPSF的收敛准则 tol 1e-7
    hasConverged = max(norm(Temp1, 'fro'), norm(Temp2, 'fro'))/M_fnorm < tol;  
    % 这个是为了和canyi比较自己设置的收敛准则 tol 1e-8
%     hasConverged = max(max(abs(Temp1(:))),max(abc(Temp2(:)))) < tol; 
end

L = X*H*Y';

function new_H = get_H_update
%Calculates a new value of H

new_H = X'*svdthresh((M - S + W + Z/mu + X*(E - N/mu)*Y')/2, 0.5/mu)*Y;
end

function new_E = get_E_update
%Calculates a new value of E

new_E = svdthresh(H - X'*W*Y + N/mu, kappa/mu);
end

function new_S = get_S_update
%Calculates a new value of S

new_S = softhresh(M - X*H*Y' + Z/mu, lambda/mu);
end

function new_Z = get_Z_update(Dir)
%Calculates a new value of Z

new_Z = Z + mu * Dir;
end

function new_N = get_N_update(Dir)
%Calculates a new value of N

new_N = N + mu * Dir;
end

end