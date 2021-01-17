function [L, S] = tpcps(M, W, varargin)

[n1, n2, n3] = size(M);
M_fnorm = norm(M(:),'fro');

if nargin > 3
    lambda = varargin{1};
    kappa = varargin{2};
else
    lambda = 1/sqrt(n3*max(n1,n2));
    if nargin > 2
        kappa = varargin{1};
    else
        kappa = 0.5 ;  
    end
end
    
mu = 1e-4; 
alpha = 1.1;
mu_max = 1e18;

L = zeros(n1, n2, n3);
E = L;
Z = L;
N = L;

tol = 1e-7;   

iter = 0;
iter_max = 1000;

hasConverged = false;

while ~hasConverged && iter < iter_max
    iter = iter + 1;
    
    S = get_S_update;
    L = get_L_update;
    E = get_E_update;
    
    Temp1 = M - S - L;
    Temp2 = L - E - W;
    
    Z = get_Z_update(Temp1);
    N = get_N_update(Temp2);
    
    mu = min(mu*alpha, mu_max);  
    
    hasConverged = max(norm(Temp1(:), 'fro'), norm(Temp2(:), 'fro'))/M_fnorm < tol; 
end

function new_L = get_L_update

[new_L,~] = prox_tnn((M - S + W + E + (Z - N)/mu)/2, 0.5/mu);

end

function new_E = get_E_update

[new_E,~] = prox_tnn(L - W + N/mu, kappa/mu);

end

function new_S = get_S_update

new_S = prox_l1(M - L + Z/mu, lambda/mu);

end

function new_Z = get_Z_update(Dir)

new_Z = Z + mu*Dir;

end

function new_N = get_N_update(Dir)

new_N = N + mu*Dir;

end

end