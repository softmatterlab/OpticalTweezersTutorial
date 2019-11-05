% PSD analysis for Nexp

% Initialization of the workspace
clear;

close all;

addpath psd



% Load data file



load('Data_positions_Fig9_1P6_S.mat')


gamma=6*pi*eta*a;

kb=1.38064852e-23;

D0=kb*T/gamma;

subs=1; %use a subsampled data set

nw=round(size(x(1:subs:end,:),1)/500); % number of windows
%[fc_psd,D_psd,sigma_fc_psd,sigma_D_psd,f,XX,fw_mean,Pk,fcut,h]=psdfit_analytic(x(1:subs:end,:),dt*subs,nw,1/4);
[fc_psd,D_psd,sigma_fc_psd,sigma_D_psd,f,XX,fw_mean,Pk,sigma_Pk,fcut]=psd_lfit(x(1:subs:end,:),dt*subs,nw,1/4);


gamma_psd=kb*T./mean(D_psd);

sigma_gamma_psd=kb*T./D_psd^2*sigma_D_psd;


sigma_k_psd=2*pi*gamma*sigma_fc_psd;

% estimation of k using the estimated gamma
k_psd=2*pi*gamma_psd.*fc_psd;

sigma_k_psd=2*pi*gamma_psd.*sigma_fc_psd+2*pi*fc_psd*sigma_gamma_psd;

disp('................')
disp('PSD analysis linear fit')

disp(['D: ' num2str(D_psd*1e12) '+-' num2str(sigma_D_psd*1e12) ' m^2/s'])

disp(['gamma: ' num2str(gamma_psd*1e9) '+-' num2str(sigma_gamma_psd*1e9) ' Ns/m'])

disp(['k: ' num2str(k_psd*1e6) '+-' num2str(sigma_k_psd*1e6) ' N/m'])


