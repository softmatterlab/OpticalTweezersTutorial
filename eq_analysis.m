% Equipartition analysis for Nexp

% Initialization of the workspace
clear;

close all;

addpath eq
load('Data_positions_Fig9_1P6_S.mat');
  
subs=3; %use a subsampled data set
[k_eq,sigma_k_eq]=eq1d(x(1:subs:1000000,:),T);

disp('................')
disp('Equipartition analysis')

disp(['k_eq: ' num2str(k_eq*1e6) '+-' num2str(sigma_k_eq*1e6) ' N/m'])




