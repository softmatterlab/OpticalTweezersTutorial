% Auto correlation analysis for Nexp

% Initialization of the workspace
close all;

clear all;

%close all;

addpath acf
addpath wlsice

    
load('Data_positions_Fig9_1P2_S.mat');



kb=1.38064852e-23;

subs=3; %use a subsampled data set
%non linear fit

[k_acf, sigma_k_acf, D_acf, sigma_D_acf, tau, mc, Ec, indc]=acf_lfit(x(1:subs:size(x,1),:),T,dt*subs);
disp('Autocorrelation function analysis by linear fitting')
 
disp(['k_acf: ' num2str(k_acf*1e6) '+-' num2str(sigma_k_acf*1e6) ' pN/um'])
 
disp(['D_acf: ' num2str(D_acf*1e12) '+-' num2str(sigma_D_acf*1e12) ' um^2/s'])
 
%disp(['gamma_acf:' num2str(gamma_acf*1e9) '+-'  num2str(sigma_gamma_acf*1e9) ' pN s/um ']);
%disp(['tau_0:' num2str(tau0_exp_lf_I*1e3) ' ms']);
