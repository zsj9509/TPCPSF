function [L, S] = tpcpsf(M, W, X, Y, varargin)

[n1, n2, n3] = size(M);
d = [size(X,2),size(Y,2)];
M_fnorm = norm(M(:),'fro');

if nargin > 5
    lambda = varargin{1};
    kappa = varargin{2};
else
    lambda = 1/sqrt(n3*max(n1,n2));
    if nargin > 4
        kappa = varargin{1};
    else
        kappa = 0.5; 
    end
end

mu = 1e-3;
alpha = 1.1;
mu_max = 1e18;

Z = zeros(n1, n2, n3);
S = Z;

H = zeros([d,n3]);
E = H;
N = H;

tol=1e-7;

iter = 0;
iter_max = 1000;

hasConverged = false;

while ~hasConverged && iter < iter_max
    iter = iter + 1;
    
    S = get_S_update;
    H = get_H_update;
    E = get_E_update;
    
    Temp1 = M - S - tprod(tprod(X,H),tran(Y));
    Temp2 = H - E - tprod(tprod(tran(X),W),Y);
    
    Z = get_Z_update(Temp1);
    N = get_N_update(Temp2);
    
    mu = min(mu*alpha, mu_max);  
    
    hasConverged = max(norm(Temp1(:), 'fro'), norm(Temp2(:), 'fro'))/M_fnorm < tol;
end

L = tprod(tprod(X,H),tran(Y));

function new_H = get_H_update
    
[new_H,~] = prox_tnn((M - S + W + Z/mu + tprod(tprod(X,(E - N/mu)),tran(Y)))/2, 0.5/mu);
new_H = tprod(tprod(tran(X),new_H),Y);

end

function new_E = get_E_update

[new_E,~] = prox_tnn(H - tprod(tprod(tran(X),W),Y) + N/mu, kappa/mu);

end

function new_S = get_S_update

new_S = prox_l1(M - tprod(tprod(X,H),tran(Y)) + Z/mu,lambda/mu);

end

function new_Z = get_Z_update(Dir)

new_Z = Z + mu * Dir;

end

function new_N = get_N_update(Dir)

new_N = N + mu * Dir;

end

end