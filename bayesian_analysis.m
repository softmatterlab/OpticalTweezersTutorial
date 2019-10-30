
clear all 
close all
%load data files
load('Data_positions_Fig9_1P2_S.mat')
addpath bayesian
xx = reshape(x,[size(x,1)*size(x,2),1]);
N=length(xx);
%some useful calculations from the positions
subs=10;

[k, sigma_k, gamma, sigma_gamma, D, sigma_D]= bayesian(xx(1:subs:end), dt*subs,T, a);
