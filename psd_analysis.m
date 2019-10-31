% PSD analysis for Nexp

% Initialization of the workspace
clear;

close all;

addpath psd

nfile='figs_psd';

% Load data file



 


load('Data_positions_Fig9_1P2_S.mat')


gamma=6*pi*eta*a;

kb=1.38064852e-23;

D0=kb*T/gamma;

subs=1; %use a subsampled data set

nw=round(size(x(1:subs:end,:),1)/500); % number of windows

[mfc_psd,mD_psd,Efc_psd,ED_psd,f,XX,fw_mean,Pk,EPk,fcut]=psd_lfit(x(1:subs:1e4,:),dt*subs,nw,1/4);

%[mfc_psd,mD_psd,Efc_psd,ED_psd,f,XX,fw_mean,Pk,EPk,fcut]=psd_lfit(Vx(1:subs:end,:),dt*subs,nw,1/4);


mgamma_psd=kb*T./mD_psd;

Egamma_psd=kb*T./mD_psd^2*ED_psd;

mk_psd=2*pi*gamma*mfc_psd;

Ek_psd=2*pi*gamma*Efc_psd;

% estimation of k using the estimated gamma
mk2_psd=2*pi*mgamma_psd.*mfc_psd;

Ek2_psd=2*pi*mgamma_psd.*Efc_psd+2*pi*mfc_psd*Egamma_psd;

disp('...')

disp('PSD analysis')

disp(['D: ' num2str(mD_psd) '+-' num2str(ED_psd) ' m^2/s'])

disp(['D0: ' num2str(D0) ' m^2/s']);

disp(['gamma: ' num2str(mgamma_psd) '+-' num2str(Egamma_psd) ' Ns/m'])

disp(['gamma0: ' num2str(gamma) ' Ns/m'])

disp(['k: ' num2str(mk_psd) '+-' num2str(Ek_psd) ' N/m'])

disp(['k2: ' num2str(mk2_psd) '+-' num2str(Ek2_psd) ' N/m'])

disp(['fc: ' num2str(mfc_psd) ' Hz'])


