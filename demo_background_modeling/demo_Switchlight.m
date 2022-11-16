addpath(genpath(cd))
clear
close all
clc
load('switchlight.mat');
sp=[25 27 40 45 168 216 226 237 246 249 277 287 289 294];
srcDir = uigetdir('.\groundtruth\groundtruthswitchlight');
cd(srcDir);
namelist=dir('*.bmp');
cd ..
len=14;
s=cell(len,1);
for i=1:len
    s{i}=imbinarize(rgb2gray(double(imread(namelist(i).name))));
end
T = double(image_M)/255;
sizT = size(T);
T_4 = T;
num = 249;
figure,imshow(T_4(:,:,:,num),'border','tight','initialmagnification','fit');
axis normal
figure,imshow(s{10},'border','tight','initialmagnification','fit');
axis normal
clear srcDir namelist I i T

F_tpcps_list=zeros(len,1);
T_W = repmat(mean(T_4(:,:,:,[1:14,65:110]),4),[1,1,1,sizT(4)]);
%%TPCPS(TR)
kappa = 2; 
[L_tpcps_tr, S_tpcps_tr] = tpcps_p_order(T_4, T_W, kappa);
for i=1:len
    J=imbinarize(rgb2gray(abs(S_tpcps_tr(:,:,:,sp(i)))));
    [F_tpcps,~] = findFMeasure(J, s{i});
    F_tpcps_list(i,1)=F_tpcps;
end
clear J F_tpcps Ori S0rgb bg L_tpcps_bg
mean(F_tpcps_list)
figure,imshow(L_tpcps_tr(:,:,:,num),'border','tight','initialmagnification','fit');
axis normal
S_tpcps_picshow=imbinarize(rgb2gray(abs(S_tpcps_tr(:,:,:,num))));
figure,imshow(S_tpcps_picshow,'border','tight','initialmagnification','fit');
axis normal

%%TPCPS(MR)
kappa = 2; 
[L_tpcps_mr, S_tpcps_mr] = tpcps_p_order_mtnn(T_4, T_W, kappa);
for i=1:len
    J=imbinarize(rgb2gray(abs(S_tpcps_mr(:,:,:,sp(i)))));
    [F_tpcps,~] = findFMeasure(J, s{i});
    F_tpcps_list(i,1)=F_tpcps;
end
clear J F_tpcps Ori S0rgb bg L_tpcps_bg
mean(F_tpcps_list)
figure,imshow(L_tpcps_mr(:,:,:,num),'border','tight','initialmagnification','fit');
axis normal
S_tpcps_picshow=imbinarize(rgb2gray(abs(S_tpcps_mr(:,:,:,num))));
figure,imshow(S_tpcps_picshow,'border','tight','initialmagnification','fit');
axis normal