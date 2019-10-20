% Auto correlation analysis for Nexp

% Initialization of the workspace
close all

clear;

%close all;

addpath msd



subs=1; %use a subsampled data set


    
load(['Data_positions_Fig9_1P2_S.mat']);


    
[k_msd, Ek_msd, D_msd, ED_msd, tau, mmsd, Emsd,indc]=msd_nfilt_covmat(x(1:subs:size(x,1),1:2),T,dt*subs,1000);


disp('')


disp(['tau0_msd: ' num2str(tau) ])%'+-' num2str(Etau)])

disp(['k_msd: ' num2str(k_msd*1e6) '+-' num2str(Ek_msd*1e6)])

disp(['D_msd: ' num2str(D_msd*1e12) '+-' num2str(ED_msd*1e12)])

%disp(['gamma_msd: ' num2str(gamma_msd*1e6) '+-' num2str(Egamma_msd*1e6)])