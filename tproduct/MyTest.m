clear;clc
% % verify trace(\mathcal{A}) = trace(\mathcal{A}^T)
% n1 = 4;
% n2 = 5;
% n3 = 6;
% r = 2;
% A = rand(n1,n1,n3);
% At= tran(A);
% trace(A(:,:,1))== trace(At(:,:,1)) % satisfy the first suppose

% % verity ||A||_F^2 = trace(A*A^T(:,:,1)) = trace(A^T*A(:,:,1))
% n1 = 4;
% n2 = 5;
% n3 = 6;
% r = 2;
% A = rand(n1,r,n3);
% At= tran(A);
% Af1 = tprod(A,At);
% Af2 = tprod(At,A);
% norm(A(:),'fro') == sqrt(trace(Af1(:,:,1)))
% trace(Af1(:,:,1))== trace(Af2(:,:,1)) % satisfy the econd supposes

% verity trace(A*B*C(:,:,1)) = trace(C*B*A(:,:,1))
n1 = 4;
n2 = 5;
n3 = 6;
r = 2;
A = rand(n1,r,n3);
B = rand(r,n2,n3);
C = rand(n2,n1,n3);
prod_ABC = tprod(tprod(A,B),C);
prod_CAB = tprod(tprod(C,A),B);
trace(prod_ABC (:,:,1))== trace(prod_CAB(:,:,1)) % satisfy the econd supposes