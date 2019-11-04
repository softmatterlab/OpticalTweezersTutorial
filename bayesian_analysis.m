
clear all 
close all
%load data files
load('Data_positions_Fig9_1P2_S.mat')
addpath bayesian
xx = reshape(x,[size(x,1)*size(x,2),1]);
N=length(xx);
%some useful calculations from the positions
subs=10;

[k_bay, sigma_k_bay, gamma_bay, sigma_gamma_bay, D_bay, sigma_D_bay]= bayesian(xx(1:subs:end), dt*subs,T, a);
disp('................')

disp('Bayesian Inference analysis')
 
disp(['k_acf: ' num2str(k_bay*1e6) '+-' num2str(sigma_k_bay*1e6) ' pN/um'])
 
disp(['D_acf: ' num2str(D_bay*1e12) '+-' num2str(sigma_D_bay*1e12) ' um^2/s'])
 
disp(['gamma_acf:' num2str(gamma_bay*1e9) '+-'  num2str(sigma_gamma_bay*1e9) ' pN s/um ']);
 
disp('................')