% Potential analysis for Nexp

%% Initialization of the workspace
clear;

close all;

addpath ../data/

% Load data file

load('Data_x_positions_Exp_III.mat');


%Boltzmann constant
kb=1.38064852e-23;

%number of bins of the histogram, if not set default is 50
P=50; 

%use a subsampled data set
subs=1;

%linear fit
disp('................')
disp('Potential analysis linear fit')
[k_pot_lf, sigma2_k_pot_lf, x_alpha_lf, mrho_lf, sigma2_rho_lf, mU_lf, sigma2_Ulf, rho0_lf, x_eq_lf, U_0_lf]=pot_lfit(x(1:subs:size(x,1),:),T,P);
disp('lineal')
disp(['k_pot: ' num2str(k_pot_lf*1e6) '+-' num2str(sigma2_k_pot_lf*1e6) ' pN/um']);
%non-linear fit
[k_potnnl, sigma2_k_potmnl, x_alpha_nf, mrho_nl, sigma2_rho_nl, mU_nl, sigma2_U_nl, rho0_nl, x_eq_nl, coefa]=pot_nlfit(x(1:subs:size(x,1),:),T,P);
disp('no lineal')
disp(['k_pot: ' num2str(k_potnnl*1e6) '+-' num2str(sigma2_k_potmnl*1e6) ' pN/um']);
