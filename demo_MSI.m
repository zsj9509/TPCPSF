addpath(genpath(cd))
clear
close all
clc

%tomatoes
cur_data='fake_and_real_tomatoes.mat';
%%TPCPSF(TR)
kappa=0.7;
atom=500;
MPSNR=MSI_denoising_tr(cur_data,kappa,atom);

%%TPCPSF(MR)
kappa=0.7;
atom=500;
MPSNR=MSI_denoising_mr(cur_data,kappa,atom);