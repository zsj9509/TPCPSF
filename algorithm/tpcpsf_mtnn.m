function [L, S] = tpcpsf_mtnn(M, W, X, Y, varargin)

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

mu = (1e-3)*n3;
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
lambda=lambda*n3;
while ~hasConverged && iter < iter_max
    iter = iter + 1;
    
    S = get_S_update;
    s=S(:,:,1);
    H = get_H_update;
    h=H(:,:,1);
    E = get_E_update;
    e=E(:,:,1);
    Temp1 = M - S - tprod(tprod(X,H),tran(Y));
    temp1=Temp1(:,:,1);
    Temp2 = H - E - tprod(tprod(tran(X),W),Y);
    temp2=Temp2(:,:,1);
    Z = get_Z_update(Temp1);
    z=Z(:,:,1);
    N = get_N_update(Temp2);
    n=N(:,:,1);
    
    mu = min(mu*alpha, mu_max);  
    
    hasConverged = max(norm(Temp1(:), 'fro'), norm(Temp2(:), 'fro'))/M_fnorm < tol;
end

L = tprod(tprod(X,H),tran(Y));
l=L(:,:,1);
function new_H = get_H_update
    
[new_H,~] = prox_tnn((M - S + W + Z/mu + tprod(tprod(X,(E - N/mu)),tran(Y)))/2, 0.5*n3/mu);
h=new_H(:,:,1);
new_H = tprod(tprod(tran(X),new_H),Y);
h1=new_H(:,:,1);

end

function new_E = get_E_update
a=tran(X);a1=a(:,:,1);b=tprod(tran(X),W);b1=b(:,:,1);
[new_E,~] = prox_tnn(H - tprod(tprod(tran(X),W),Y) + N/mu, n3*kappa/mu);
e=new_E(:,:,1);
end

function new_S = get_S_update
a=tran(X);a1=a(:,:,1);b=tprod(tran(X),W);b1=b(:,:,1);
new_S = prox_l1(M - tprod(tprod(X,H),tran(Y)) + Z/mu,lambda/mu);
s=new_S(:,:,1);
end

function new_Z = get_Z_update(Dir)

new_Z = Z + mu * Dir;

end

function new_N = get_N_update(Dir)

new_N = N + mu * Dir;

end

end