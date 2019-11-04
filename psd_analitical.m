% PSD analysis for Nexp

% Initialization of the workspace
clear all;

close all;

addpath psd



% Load data file



load('Data_positions_Fig9_1P6_S.mat')


gamma=6*pi*eta*a;

kb=1.38064852e-23;

D0=kb*T/gamma;

subs=3; %use a subsampled data set

nw=round(size(x(1:subs:end,:),1)/500);
%%ANALITICAL FIT

[fc_psd,D_psd,sigma_fc_psd,sigma_D_psd,f,XX,fw_mean,Pk,fcut,h]=psdfit_analytic(x(1:subs:end,:),dt*subs,nw,1/4);
%%




gamma_psd=kb*T./mean(D_psd);
sigma_D_var=std(D_psd)
sigma_gamma_psd=kb*T./mean(D_psd)^2*sigma_D_var;

% estimation of k using the estimated gamma
k_psd=2*pi*gamma_psd.*mean(fc_psd);
sigma_fc_var=std(fc_psd);
sigma_k_psd=2*pi*(gamma_psd.*sigma_fc_var+mean(fc_psd)*sigma_gamma_psd);

disp('................')
disp('PSD analitycal solution')

disp(['D: ' num2str(mean(D_psd)*1e12) '+-' num2str(sigma_D_var*1e12) ' um^2/s'])

disp(['gamma: ' num2str(gamma_psd*1e9) '+-' num2str(sigma_gamma_psd*1e9) ' pNms/um'])

disp(['k: ' num2str(k_psd*1e6) '+-' num2str(sigma_k_psd*1e6) ' pN/um'])



