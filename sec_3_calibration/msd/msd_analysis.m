% Auto correlation analysis for Nexp

% Initialization of the workspace
clear all

clear;

%close all;

addpath ../data/
addpath ../statistics_func/

subs=1; %use a subsampled data set

maxlag=300;
    
load(['Data_x_positions_Exp_I.mat']);

kB=1.38064852e-23;
[k_msd,sigma_k_msd, tau0, sigma2_tau0, D_msd, sigma_D_msd, tau, mmsd, sigma_msd, indc, gamma_msd, sigma2_gamma_msd] =msd_nlfit(x(1:subs:size(x,1),:),T, dt*subs,maxlag);   

gamma_msd=kB*T/D_msd;
sigma_gamma_msd=kB*T/D_msd^2*sigma_D_msd;
disp('................')


disp(['tau0_msd: ' num2str(tau0) ])%'+-' num2str(Etau)])

disp(['k_msd: ' num2str(k_msd*1e6) '+-' num2str(sigma_k_msd*1e6)])

disp(['D_msd: ' num2str(D_msd*1e12) '+-' num2str(sigma_D_msd*1e12)])


disp(['gamma_msd: ' num2str(gamma_msd*1e9) '+-' num2str(sigma_gamma_msd*1e9)])


