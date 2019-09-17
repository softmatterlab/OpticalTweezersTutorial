% Equipartition analysis for Nexp

% Initialization of the workspace
clear;

%close all;

addpath eq
load(['Data_positions_Fig9_1P4_S.mat']);
  
subs=1; %use a subsampled data set
[k_eq,Ek_eq]=eq1d_v1(x(1:subs:1000000,:),T,0e-9);



