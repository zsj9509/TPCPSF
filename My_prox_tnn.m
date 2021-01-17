function [X,tnn,trank] = My_prox_tnn(Y,rho)

sa = size(Y);
la = length(sa);
tnn = 0;
trank = 0;

for i = 3:la
    Y = fft(Y,[],i);
end

faces = prod(sa(3:la));
sU = sa;         
sU(2) = sU(1);
sV = sa;        
sV(1) = sV(2);
U = zeros(sU); 
S = zeros(sa);
V = zeros(sV);
[U,S,V] = takeSVDs(U,S,V,Y,faces);
S_true = zeros(sa);
X = zeros(sa);

for i = 1:faces
    dS = S(:,:,i);
    S_diag = SoftThresh(diag(dS),rho);
    r = length(S_diag);
    if r>=1
        S_true(1:r,1:r,i) = diag(S_diag);
        X(:,:,i) = U(:,1:r,i)*diag(S_diag)*V(:,1:r,i)'; 
        tnn = tnn+sum(S_diag);
        trank = max(trank,r);
    end
end
for i = la:-1:3
    X = ifft(X,[],i);
end
X = real(X);
tnn = sum(diag( S_true(:,:,1) ) );
end