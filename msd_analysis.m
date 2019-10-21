% Auto correlation analysis for Nexp

% Initialization of the workspace
close all

clear;

%close all;

addpath msd



subs=10; %use a subsampled data set
%exp I SUBS=20, maxlag=50 
%exp II subs-17, maxlag=25
%exp III subs=10, maxlag=25
maxlag=25;
    
load(['Data_positions_Fig9_1P6_S.mat']);


    
[k_msd, Ek_msd, D_msd, ED_msd, tau, mmsd, Emsd,indc]=msd_nfilt(x(1:subs:size(x,1),:),T,dt*subs,maxlag);


disp('')


disp(['tau0_msd: ' num2str(tau) ])%'+-' num2str(Etau)])

disp(['k_msd: ' num2str(k_msd*1e6) '+-' num2str(Ek_msd*1e6)])

disp(['D_msd: ' num2str(D_msd*1e12) '+-' num2str(ED_msd*1e12)])

%disp(['gamma_msd: ' num2str(gamma_msd*1e6) '+-' num2str(Egamma_msd*1e6)])