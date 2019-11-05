
clear all 
close all
%load data files
load('Data_positions_Fig9_1P6_S.mat')
x = x - repmat(mean(x),size(x,1),1);
addpath bayesian
xx = reshape(x,[size(x,1)*size(x,2),1]);
N=length(xx);
%some useful calculations from the positions
nsubs=3;
subs=1;
[k_bay, sigma_k_bay, gamma_bay, sigma_gamma_bay, D_bay, sigma_D_bay]= bayesian(xx(1:subs:end),dt*subs,T, a, nsubs);
%[k_bay, sigma_k_bay, gamma_bay, sigma_gamma_bay, D_bay, sigma_D_bay]= bayesian([],dt,T, a);
disp('................')

disp('Bayesian Inference analysis')
 
disp(['k_bay: ' num2str(k_bay*1e6) '+-' num2str(sigma_k_bay*1e6) ' pN/um'])
 
disp(['D_bay: ' num2str(D_bay*1e12) '+-' num2str(sigma_D_bay*1e12) ' um^2/s'])
 
disp(['gamma_bay:' num2str(gamma_bay*1e9) '+-'  num2str(sigma_gamma_bay*1e9) ' pN ms/um ']);
 
disp('................')