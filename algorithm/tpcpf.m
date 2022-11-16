function [L, S] = tpcpf(M, X, Y, varargin)

[n1, n2, n3] = size(M);
d = [size(X,2),size(Y,2)];
M_fnorm = norm(M(:),'fro');

if nargin > 4
    lambda = varargin{1};
else
    lambda = 1/sqrt(n3*max(n1,n2));
end

mu = 1e-3;   
alpha = 1.1;
mu_max = 1e18;

Z = zeros(n1, n2, n3);
S = Z; 
H = zeros([d,n3]);

tol = 1e-7;   

iter = 0;
iter_max = 1000;

hasConverged = false;

while ~hasConverged && iter < iter_max
    iter = iter + 1;
    
    S = get_S_update;
    H = get_H_update;
    
    Temp1 = M - S - tprod(tprod(X,H),tran(Y));
    
    Z = get_Z_update(Temp1);
    
    mu = min(mu*alpha, mu_max);  
    
    hasConverged = norm(Temp1(:), 'fro')/M_fnorm < tol; 
end

L = tprod(tprod(X,H),tran(Y));

function new_H = get_H_update

[new_H,~] = prox_tnn((M - S + Z/mu + tprod(tprod(X,H),tran(Y)))/2, 0.5/mu);
new_H = tprod(tprod(tran(X),new_H),Y);

end

function new_S = get_S_update
    
new_S = prox_l1(M - tprod(tprod(X,H),tran(Y)) + Z/mu,lambda/mu);

end

function new_Z = get_Z_update(Dir)

new_Z = Z + mu * Dir;

end

end