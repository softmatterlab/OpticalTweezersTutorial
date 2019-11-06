% Auto correlation analysis for Nexp

% Initialization of the workspace
close all;

clear all;

%close all;

addpath acf
addpath wlsice

    
load('Data_positions_Fig9_1P6_S.mat');



kB=1.38064852e-23;

subs=1; %use a subsampled data set
%non linear fit



[k_acf_lf, sigma_k_acf_lf, D_acf_lf, sigma_D_acf_lf,gamma_acf_lf, sigma_gamma_acf_lf,tau_acf_lf, mc, Ec,indc, tau0_exp_lf, c0_exp_lf]=acf_lfit(x(1:subs:size(x,1),:),T,dt*subs);

[k_acf_nl, sigma_k_acf_nl, D_acf_nl, sigma_D_acf_nl,gamma_acf_nl, sigma_gamma_acf_nl, tau_nl, mc, Ec, indc, tau0_exp_nl, c0_exp_nl]=acf_nlfit(x(1:subs:size(x,1),:),T,dt*subs);

disp('................')
disp('Autocorrelation function analysis by linear fitting')
 
disp(['k_acf: ' num2str(k_acf_lf*1e6) '+-' num2str(sigma_k_acf_lf*1e6) ' pN/um'])
 
disp(['D_acf: ' num2str(D_acf_lf*1e12) '+-' num2str(sigma_D_acf_lf*1e12) ' um^2/s'])

disp(['gamma_acf:' num2str(gamma_acf_lf*1e9) '+-'  num2str(sigma_gamma_acf_lf*1e9) ' pN ms/um ']);



disp('................')
disp('Autocorrelation function analysis by non-linear fitting')
 
disp(['k_acf: ' num2str(k_acf_nl*1e6) '+-' num2str(sigma_k_acf_nl*1e6) ' pN/um'])
 
disp(['D_acf: ' num2str(D_acf_nl*1e12) '+-' num2str(sigma_D_acf_nl*1e12) ' um^2/s'])

disp(['gamma_acf:' num2str(gamma_acf_nl*1e9) '+-'  num2str(sigma_gamma_acf_nl*1e9) ' pN ms/um ']);
