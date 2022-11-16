addpath(genpath(cd))
clear
close all
clc
load('Curtain.mat');
sp=[16 18 59];
srcDir = uigetdir('.\groundtruth\groundtruthcurtain');
cd(srcDir);
namelist=dir('*.bmp');
cd ..
len=3;
s=cell(len,1);
for i=1:len
    s{i}=imbinarize(rgb2gray(double(imread(namelist(i+2).name))));
end
T = double(I(:,:,:,151:210))/255;
maxT = max(abs(T(:)));
sizT = size(T);
T_4 = T;
n1 = sizT(1)*sizT(2);
n2 = sizT(3);
n3 = sizT(4);
T_3 = reshape(T,[n1,n2,n3]);
num = 18;
figure,imshow(T_4(:,:,:,num),'border','tight','initialmagnification','fit');
axis normal
figure,imshow(s{2},'border','tight','initialmagnification','fit');
axis normal
clear srcDir namelist i T

% % tpcps
% F_tpcps_list=zeros(len,1);
% tmp= double(I)/255;
% T_W = repmat(mean(tmp(:,:,:,1:73),4),[1,1,1,60]);
% %% tpcps(tr)
% kappa = 10;
% [L_tpcps_tr, S_tpcps_tr] = tpcps_p_order(T_4, T_W, kappa);
% for i=1:len
%     J=imbinarize(rgb2gray(abs(S_tpcps_tr(:,:,:,sp(i)))));
%     [F_tpcps,~] = findFMeasure(J, s{i});
%     F_tpcps_list(i,1)=F_tpcps;
% end
% clear J F_tpcps Ori S0rgb bg L_tpcps_bg
% mean(F_tpcps_list)
% figure,imshow(L_tpcps_tr(:,:,:,num),'border','tight','initialmagnification','fit');
% axis normal
% S_tpcps_picshow=imbinarize(rgb2gray(abs(S_tpcps_tr(:,:,:,num))));
% figure,imshow(S_tpcps_picshow,'border','tight','initialmagnification','fit');
% axis normal
% %%tpcps(mr)
% kappa = 10;
% [L_tpcps_mr, S_tpcps_mr] = tpcps_p_order_mtnn(T_4, T_W, kappa);
% for i=1:len
%     J=imbinarize(rgb2gray(abs(S_tpcps_mr(:,:,:,sp(i)))));
%     [F_tpcps,~] = findFMeasure(J, s{i});
%     F_tpcps_list(i,1)=F_tpcps;
% end
% clear J F_tpcps Ori S0rgb bg L_tpcps_bg
% mean(F_tpcps_list)
% figure,imshow(L_tpcps_mr(:,:,:,num),'border','tight','initialmagnification','fit');
% axis normal
% S_tpcps_picshow=imbinarize(rgb2gray(abs(S_tpcps_mr(:,:,:,num))));
% figure,imshow(S_tpcps_picshow,'border','tight','initialmagnification','fit');
% axis normal

%% tpcp
opts.tol = 1e-7; 
opts.mu = 1e-4;
opts.rho = 1.1;
opts.DEBUG = 1;
opts.max_iter = 1000;
lambda = 1/sqrt(max(n1,n2)*n3);
[Xhat,E,~,~] = trpca_tnn(T_3,lambda,opts);%60 
Xhat = max(Xhat,0);
Xhat = min(Xhat,maxT);
F_tpcp_list=zeros(len,1);
for i=1:len
    J=imbinarize(rgb2gray(abs(reshape(E(:,:,sp(i)),[sizT(1) sizT(2) sizT(3)]))));
    [F_tpcp,~] = findFMeasure(J, s{i});
    F_tpcp_list(i,1)=F_tpcp;
end
mean(F_tpcp_list)
clear J i opts
L_tpcp_picshow=reshape(Xhat(:,:,num),[sizT(1) sizT(2) sizT(3)]);
figure,imshow(L_tpcp_picshow,'border','tight','initialmagnification','fit');
axis normal
S_tpcp_picshow=imbinarize(rgb2gray(abs(reshape(E(:,:,num),[sizT(1) sizT(2) sizT(3)]))));
figure,imshow(S_tpcp_picshow,'border','tight','initialmagnification','fit');
axis normal

%% ÌôÆäËûµÄ60Ö¡
% addpath(genpath(cd))
% clear
% close all
% clc
% load('Curtain.mat');
% sp=[91 93 166 168 209];
% srcDir = uigetdir('.\groundtruth\groundtruthcurtain');
% cd(srcDir);
% namelist=dir('*.bmp');
% cd ..
% len=2;
% s=cell(len,1);
% for i=1:len
%     s{i}=imbinarize(rgb2gray(double(imread(namelist(i).name))));
% end
% T = double(I)/255;
% sizT = size(T);
% T_4 = T;
% T_W = repmat(mean(T_4(:,:,:,1:73),4),[1,1,1,60]);
% clear srcDir namelist I i T
% kappa = 9;
% T_4_60=T_4(:,:,:,83:142); %150:209
% num=9;
% figure,imshow(T_4_60(:,:,:,num),'border','tight','initialmagnification','fit');
% axis normal
% figure,imshow(s{1},'border','tight','initialmagnification','fit');
% axis normal
% [L_tpcps, S_tpcps] = tpcps_p_order(T_4_60, T_W, kappa);
% F_tpcps_list=zeros(len,1);
% for i=1:len
%     J=imbinarize(rgb2gray(abs(reshape(S_tpcps(:,:,:,sp(i)-83+1),[sizT(1) sizT(2) sizT(3)]))));
%     [F_tpcps,~] = findFMeasure(J, s{i});
%     F_tpcps_list(i,1)=F_tpcps;
% end
% mean(F_tpcps_list)
% L_tpcps_picshow=reshape(L_tpcps(:,:,:,num),[sizT(1) sizT(2) sizT(3)]);
% figure,imshow(L_tpcps_picshow,'border','tight','initialmagnification','fit');
% axis normal
% S_tpcps_picshow=imbinarize(rgb2gray(abs(reshape(S_tpcps(:,:,:,num),[sizT(1) sizT(2) sizT(3)]))));
% figure,imshow(S_tpcps_picshow,'border','tight','initialmagnification','fit');
% axis normal