% Auto correlation analysis for Nexp

% Initialization of the workspace
clear all

clear;

%close all;

addpath msd

addpath wlsice

subs=3; %use a subsampled data set
%exp I SUBS=20, maxlag=50 
%exp II subs-17, maxlag=25
%exp III subs=10, maxlag=25
maxlag=300;
    
load(['Data_positions_Fig9_1P6_S.mat']);

kB=1.38064852e-23;
[k_msd,sigma_k_msd, tau0, sigma2_tau0, D_msd, sigma_D_msd, tau, mmsd, sigma_msd, indc, gamma_msd, sigma2_gamma_msd] =msd_nfilt(x(1:subs:size(x,1),:),T, dt*subs,maxlag);   

gamma_msd=kB*T/D_msd;
sigma_gamma_msd=kB*T/D_msd^2*sigma_D_msd;

disp('')


disp(['tau0_msd: ' num2str(tau0) ])%'+-' num2str(Etau)])

disp(['k_msd: ' num2str(k_msd*1e6) '+-' num2str(sigma_k_msd*1e6)])

disp(['D_msd: ' num2str(D_msd*1e12) '+-' num2str(sigma_D_msd*1e12)])


disp(['gamma_msd: ' num2str(gamma_msd*1e9) '+-' num2str(sigma_gamma_msd*1e9)])