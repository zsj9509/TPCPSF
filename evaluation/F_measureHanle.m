function [Fmeasure]=F_measureHanle
%本程序的功能是对显著性特征提取的结果计算F-measure值。
%by hanlestudy@163.com
clc
clear
imnames=dir(path_output);  
imnames2=dir(path_target);  
num=length(imnames);
belt2=0.3;
reca = zeros(num,1);
prec = zeros(num,1);
for i=1:num
    Target=imread(imnames2(i).name);%读图
    mark = Target(:,:,1);
    mark = mark/max(mark(:));%二值化Ground-truth
    Output=imread(imnames(i).name);
    thresh=2*mean(mean(mean(Output)));  %自适应阈值
    label=reshape(mark,1,256*256);
    score=reshape(Output(:,:,1),1,256*256); 
    sco_th0=(score)>thresh;
    sco_th=uint8(sco_th0);
    TP = length(find((label == 1) & (sco_th == 1)));
    FP = length(find((label == 0) & (sco_th == 1)));
    FN = length(find((label == 1) & (sco_th == 0)));
    reca(i,1) = TP/(TP+FN);
    prec(i,1) = TP/(TP+FP);
    i
end
P=mean(prec);
R=mean(reca);
Fmeasure=((1+belt2)*P*R)/(belt2*P+R)
