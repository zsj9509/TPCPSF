addpath(genpath(cd))
clear
close all
clc

%% load data 以Curtain为例 300帧
load('Curtain.mat');
sp=[91 93 166 168 209];
srcDir = uigetdir('.\groundtruth\groundtruthcurtain');
cd(srcDir);
namelist=dir('*.bmp');
cd ..
len=5;
s=cell(len,1);
for i=1:len
    s{i}=imbinarize(rgb2gray(double(imread(namelist(i).name))));
end
T = double(I)/255;
T_4=T;
maxT = max(abs(T(:)));
sizT = size(T);
n1 = sizT(1)*sizT(2);
n2 = sizT(3);
n3 = sizT(4);
T = reshape(T,[n1,n2,n3]);
M = reshape(T,[n1*n2,n3]);
num = 168;
figure,imshow(T_4(:,:,:,num),'border','tight','initialmagnification','fit');
axis normal
clear srcDir namelist  i

%% pcp
% [L_pcp,S_pcp]=pcp(M);
% F_pcp_list=zeros(len,1);
% for i=1:len
%     J=imbinarize(rgb2gray(abs(reshape(S_pcp(:,sp(i)),[sizT(1) sizT(2) sizT(3)]))));
%     [F_pcp,~] = findFMeasure(J, s{i});
%     F_pcp_list(i,1)=F_pcp;
% end
% mean(F_pcp_list)
% clear J F_pcp i
% L_pcp_picshow=reshape(L_pcp(:,num),[sizT(1) sizT(2) sizT(3)]);
% figure,imshow(L_pcp_picshow,'border','tight','initialmagnification','fit');
% axis normal
% S_pcp_picshow=imbinarize(rgb2gray(abs(reshape(S_pcp(:,num),[sizT(1) sizT(2) sizT(3)]))));
% figure,imshow(S_pcp_picshow,'border','tight','initialmagnification','fit');
% axis normal
% 
% %% pcps
% kappa=0.8;
% W = repmat(mean(M(:,1:73),2),[1,n3]);
% [L_pcps, S_pcps]=pcps(M, W, kappa);
% F_pcps_list=zeros(len,1);
% for i=1:len
%     J=imbinarize(rgb2gray(abs(reshape(S_pcps(:,sp(i)),[sizT(1) sizT(2) sizT(3)]))));
%     [F_pcps,~] = findFMeasure(J, s{i});
%     F_pcps_list(i,1)=F_pcps;
% end
% mean(F_pcps_list)
% L_pcps_picshow=reshape(L_pcps(:,num),[sizT(1) sizT(2) sizT(3)]);
% figure,imshow(L_pcps_picshow,'border','tight','initialmagnification','fit');
% axis normal
% S_pcps_picshow=imbinarize(rgb2gray(abs(reshape(S_pcps(:,num),[sizT(1) sizT(2) sizT(3)]))));
% figure,imshow(S_pcps_picshow,'border','tight','initialmagnification','fit');
% axis normal
% clear J F_pcps i

%% tpcp
% opts.tol = 1e-7; 
% opts.mu = 1e-4;
% opts.rho = 1.1;
% opts.DEBUG = 1;
% opts.max_iter = 1000;
% lambda = 1/sqrt(max(n1,n2)*n3);
% [Xhat,E,err,iter] = trpca_tnn(T,lambda,opts);
% Xhat = max(Xhat,0);
% Xhat = min(Xhat,maxT);
% F_tpcp_list=zeros(len,1);
% for i=1:len
%     J=imbinarize(rgb2gray(abs(reshape(E(:,:,sp(i)),[sizT(1) sizT(2) sizT(3)]))));
%     [F_tpcp,~] = findFMeasure(J, s{i});
%     F_tpcp_list(i,1)=F_tpcp;
% end
% mean(F_tpcp_list)
% clear J F_tpcp i
% L_tpcp_picshow=reshape(Xhat(:,:,num),[sizT(1) sizT(2) sizT(3)]);
% figure,imshow(L_tpcp_picshow,'border','tight','initialmagnification','fit');
% axis normal
% S_tpcp_picshow=imbinarize(rgb2gray(abs(reshape(E(:,:,num),[sizT(1) sizT(2) sizT(3)]))));
% figure,imshow(S_tpcp_picshow,'border','tight','initialmagnification','fit');
% axis normal

%% brtf
[imHeight,imWidth,~,nFrames]=size(I);
IV =10;
[model] = BayesRCP(T, 'init', 'rand', 'maxRank', 10, 'maxiters', 20, 'initVar', IV, 'updateHyper', 'off',...
    'tol', 1e-3, 'dimRed', 1, 'verbose', 1);
X = double(ktensor(model.Z)) * 255;
S = model.E * 255;
X_BRCPF = reshape(X, [imHeight, imWidth, 3, nFrames])./255;
S_BRCPF = reshape(S, [imHeight, imWidth, 3, nFrames])./255;
F_BRTF_list=zeros(len,1);
for i=1:len
    J=imbinarize(rgb2gray(reshape(abs(S_BRCPF(:,:,:,sp(i))),[sizT(1) sizT(2) sizT(3)])));
    [F_BRTF,~] = findFMeasure(J, s{i});
    F_BRTF_list(i,1)=F_BRTF;
end
mean(F_BRTF_list)
L_BRTF_picshow=X_BRCPF(:,:,:,num);
figure,imshow(L_BRTF_picshow,'border','tight','initialmagnification','fit');
axis normal
S_BRTF_picshow=imbinarize(rgb2gray(abs(S_BRCPF(:,:,:,num))));
figure,imshow(S_BRTF_picshow,'border','tight','initialmagnification','fit');
axis normal
clear J i F_BRTF

% %% etrpca 使用前请从当前路径里移除 algorithm和tproduct文件夹路径，因为这两个文件夹里有函数与etrpca方法的func文件夹里函数同名
% [n1, n2, n3, n4] = size(I);
% T_3 = T;
% F_etrpca_list=zeros(len,1);
% n = min(n1*n2, n3);
% w = [];
% w = [w; 1*ones(1,1)];
% w = [w; 4*ones(n-1,1)];
% p = 1;
% [D, E, obj, err, iter]  = etrpca_tnn_lp(T_3, 1/sqrt(max(n1*n2, n3)*n4), w, p);
% E_T3 = reshape(E, [n1, n2, n3, n4]);
% D_T3 = reshape(D, [n1, n2, n3, n4]);
% for i=1:len
%     J=imbinarize(rgb2gray(abs(E_T3(:,:,:,sp(i)))));
%     [F_etrpca,~] = findFMeasure(J, s{i});
%     F_etrpca_list(i,1)=F_etrpca;
% end
% mean(F_etrpca_list)
% clear J i F_etrpca
% figure,imshow(D_T3(:,:,:,num),'border','tight','initialmagnification','fit');
% axis normal
% S_etrpca_picshow=imbinarize(rgb2gray(abs(E_T3(:,:,:,num))));
% figure,imshow(S_etrpca_picshow,'border','tight','initialmagnification','fit');
% axis normal