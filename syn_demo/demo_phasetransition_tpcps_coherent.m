% phase transition experiment  -- coherent
addpath(genpath(cd))
clear
close all
clc

n1 = 100;  % n1= 200;
n2 = n1;
n3 = n1;
ds = 50;

% trpca_result_Lr = zeros(ds,ds);
tpcps_result_Lr = zeros(ds,ds);
mkdir('./syn_tpcps_result');

num1=0; %rank
num2=0; %sparse
for rank_ratio = 0.01:0.01:0.5
    num2=0;
    num1=num1+1;
    disp(num1);
    r = ceil(rank_ratio*n1);
    L1 = randn(n1,r,n3)/sqrt(n1);
    L2 = randn(r,n2,n3)/sqrt(n2);
    L = tprod(L1,L2); % low rank part
   
    % add gaussian noise
    noise = randn(n1,n2,n3)/sqrt(n1)*sqrt(r)*1e-6;
    % noise = tprod(randn(n1,1,n3),randn(1,n2,n3))/sqrt(n1)*r*1e-5;
    W = L+noise;
    for p = 0.01:0.01:0.5
        num2=num2+1;
        m = ceil(p*n1*n2*n3);
        temp = rand(n1*n2*n3,1);
        [B,I] = sort(temp);
        I = I(1:m);
        Omega = zeros(n1,n2,n3);
        Omega(I) = 1;
        
        S = Omega.*sign(L); % sparse part, S = P_Omega sgn(L)

        M = W + S;
        

        %--------------------------- tpcps-----------------
        [L_tpcps,~] = tpcps(M,W);
        Lr_tpcps = norm(L(:)-L_tpcps(:))/norm(L(:));
        tpcps_result_Lr(num1,num2) = Lr_tpcps;

    end
    save(strcat('./syn_tpcps_result/rank',num2str(num1),'.mat'),'tpcps_result_Lr');
end

save ./syn_tpcps_result/tpcps.mat tpcps_result_Lr
