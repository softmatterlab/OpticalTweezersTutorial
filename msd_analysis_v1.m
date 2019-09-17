% Auto correlation analysis for Nexp

% Initialization of the workspace
close all

clear;

%close all;

addpath msd



subs=1; %use a subsampled data set


    
load(['Data_positions_Fig9_1P4_S.mat']);



[k_msd, Ek_msd, D_msd, ED_msd, tau, mmsd, Emsd,indc]=msd_nlfit_v1(x(1:subs:size(x,1),:),T,dt*subs,1000);

