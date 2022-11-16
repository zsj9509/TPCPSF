addpath(genpath(cd))
clear
close all
clc
load('yaleB06.mat');
M = img./255;
clear img
TM = reshape(M,[192 168 64]);
[n1,n2,n3] = size(TM);
TW = repmat(mean(TM,3),[1 1 n3]);
Xn = M;
rhos = 0.2; 
ind = find(rand(n1*n2*n3,1)<rhos);
Xn(ind) = rand(length(ind),1);
TXn = reshape(Xn,[192 168 64]);
clear ind rhos
num = 25;
raw_picshow = reshape(M(:,num),[n1 n2]);
figure,imshow(raw_picshow,'border','tight','initialmagnification','fit');
set (gcf,'Position',[0,0,168,192]);
axis normal
noisy_picshow = reshape(Xn(:,num),[n1 n2]);
figure,imshow(noisy_picshow,'border','tight','initialmagnification','fit');
set (gcf,'Position',[0,0,168,192]);
axis normal

load yaleB06_selected16.mat
Side_origion = repmat(img, [1,4]) ./ 255;
clear img
Xs = reshape(Side_origion,[n1 n2 n3]);
atom = 168;
D0 = ones(size(Xs,1), atom);
D0 = repmat(D0 / sqrt(size(Xs, 3)), [1 1 size(Xs, 3)]);    
D = ktsvd(Xs, D0, 0.1, 10);
D0Y = ones(size(Xs, 2),atom);
D0Y = repmat(D0Y / sqrt(size(Xs, 3)),[1 1 size(Xs, 3)]);    
DY = ktsvd(tran(Xs), D0Y, 0.1, 10);
tX = D;
tY = DY;
X = fft(tX,[],3);
Y = fft(tY,[],3);
tXz = zeros(size(X));
tYz = zeros(size(Y));
for i = 1:size(X,3)
    a = orth(X(:,:,i));
    [m1,m2] = size(a);
    tXz(1:m1,1:m2,i) = a;
    b = orth(Y(:,:,i));
    [m3,m4] = size(b);
    tYz(1:m3,1:m4,i) = b;
end
tXnew = ifft(tXz,[],3);
tYnew = ifft(tYz,[],3);
clear m1 m2 m3 m4 a b i tXz tYz tX tY D0 D0Y Xs Side_origion

kappa = 2;
%%TPCPSF(TR)
[L_tpcpsf_tr,S_tpcpsf_tr] = tpcpsf(TXn, TW, tXnew, tYnew, kappa);
L_tpcpsf_tr_picshow = L_tpcpsf_tr(:,:,num);
figure,imshow(L_tpcpsf_tr_picshow,'border','tight','initialmagnification','fit');
set (gcf,'Position',[0,0,168,192]);
axis normal
%%TPCPSF(MR)
[L_tpcpsf_mr,S_tpcpsf_mr] = tpcpsf_mtnn(TXn, TW, tXnew, tYnew, kappa(i));
L_tpcpsf_picshow = L_tpcpsf_mr(:,:,num);
figure,imshow(L_tpcpsf_mr_picshow,'border','tight','initialmagnification','fit');
set (gcf,'Position',[0,0,168,192]);
axis normal

