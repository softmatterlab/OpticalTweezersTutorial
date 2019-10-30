% Auto correlation analysis for Nexp

% Initialization of the workspace
close all;

clear all;

%close all;

addpath acf
addpath wlsice

    
load('Data_positions_Fig9_1P2_S.mat');



kb=1.38064852e-23;

subs=1; %use a subsampled data set
%non linear fit

[k_acf, Ek_acf, D_acf, ED_acf, tau, mc, Ec, indc]=acf_nlfit(x(1:subs:size(x,1),:),T,dt*subs);


[k_acf, Ek_acf, D_acf, ED_acf, tau, mc, Ec,indc, tau0_exp, c0_exp, acf, C]=acf_lfit(x(1:subs:size(x,1),:),T,dt*subs);

