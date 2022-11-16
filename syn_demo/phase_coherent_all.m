% figure;
clear
addpath(genpath(pwd));
load('.\coherent_all\tpcp_coherent_all.mat')
load('.\coherent_all\tpcpsf_coherent_all.mat')
load('.\coherent_all\tpcpf_coherent_all.mat')
load('.\coherent_all\tpcps_coherent_all.mat')
trpca =double(flip(tpcp_result_Lr',1)<1e-3);
idx1=find(trpca==1); 
% trpca(:,:) = 0.8;
trpca(idx1) = 0.1;

tpcps =double(flip(tpcps_result_Lr',1)<1e-3);
idx2=find(tpcps==1);
% tpcpsf(:,:) = 0.8;
tpcps(idx2) = 0.1;

tpcpf =double(flip(tpcpf_result_Lr',1)<1e-3);
idx3=find(tpcpf==1); 
% tpcpf(:,:) = 0.8;
tpcpf(idx3) = 0.1;

tpcpsf =double(flip(tpcpsf_result_Lr',1)<1e-3);
idx4=find(tpcpsf==1);
% tpcpsf(:,:) = 0.8;
tpcpsf(idx4) = 0.1;

% tpcps = double(flip(tpcps_result_Lr',1)<1e-3);
% idx4=find(tpcps==1);
% % tpcpsf(:,:) = 0.8;
% tpcps(idx4) = 0.1;

all = trpca+tpcps+tpcpf+tpcpsf;
imagesc(all);
% imagesc(trpca);
% hold on;
% imshow(tpcpf);
% hold on;
% imshow(tpcpsf);
% hold on;
axis on
axis ij
set(gca,'YtickLabel',[0.4:-0.1:0])
set(gca,'XtickLabel',[0.1:0.1:0.5])
set(gca,'FontSize',15);
xlabel('rank_t/n','fontsize',15);
ylabel('\rho_s','fontsize',15);
set(gcf,'color','w')
set(gcf,'unit','normalized','position',[0.2,0.2,0.3,0.4]);
% export_fig '.\' 'trpca_phase.pdf';