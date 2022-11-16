function [L, S] = pcp(M, varargin)
%Runs the principal component pursuit using side information algorithm
%implemented using inexact augmented Lagrange multiplier method

%%Initialisation
[m, n] = size(M);
M_fnorm = norm(M,'fro');

if nargin > 3
    lambda = varargin{1};
    kappa = varargin{2};
else
    lambda = 1/sqrt(max(m,n));
    
    if nargin > 2
        kappa = varargin{1};
    else
        kappa = 0;
    end
end
    
mu = 1/lansvd(M, 1, 'L');
alpha = 1.1;
mu_max = 1e18;

L = zeros(m, n);
E = L;
Z = L;
N = L;
W = L;

tol = 1e-7;

iter = 0;
iter_max = 1000;

hasConverged = false;

%Calculation
while ~hasConverged && iter < iter_max
    iter = iter + 1;
    
    S = get_S_update;
    L = get_L_update;
    E = get_E_update;
    
    %cache variables
    Temp1 = M - S - L;
    Temp2 = L - E - W;
    
    Z = get_Z_update(Temp1);
    N = get_N_update(Temp2);
    
    mu = min(mu*alpha, mu_max);  
    
    hasConverged = max(norm(Temp1, 'fro'), norm(Temp2, 'fro'))/M_fnorm < tol;
end

function new_L = get_L_update
%Calculates a new value of L

new_L = svdthresh((M - S + W + E + (Z - N)/mu)/2, 0.5/mu);
end

function new_E = get_E_update
%Calculates a new value of E

new_E = svdthresh(L - W + N/mu, kappa/mu);
end

function new_S = get_S_update
%Calculates a new value of S

new_S = softhresh(M - L + Z/mu, lambda/mu);
end

function new_Z = get_Z_update(Dir)
%Calculates a new value of Z

new_Z = Z + mu*Dir;
end

function new_N = get_N_update(Dir)
%Calculates a new value of N

new_N = N + mu*Dir;
end

end