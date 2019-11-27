% Equipartition analysis for Nexp

% Initialization of the workspace
clear;

close all;

addpath eq
load('Data_positions_Fig9_1P6_S.mat');
  deltax=0;
subs=1; %use a subsampled data set
[k_eq,sigma_k_eq]=eq1d(x(1:subs:size(x,1),:),T,deltax);

disp('................')
disp('Equipartition analysis considering deltax equal to zero')

disp(['k_eq: ' num2str(k_eq*1e6) '+-' num2str(sigma_k_eq*1e6) ' N/m'])

%

equipdx_tex=fopen('equipdx.txt', 'w')
load('Data_positions_Fig9_1P2_S.mat'); 

[k_eq_delta_I,sigma_k_eq_delta_I]=eq1d(x(1:subs:size(x,1),:),T,1e-9);

disp('................')
disp('Equipartition analysis considering deltax=1nm')

disp(['k_eq: ' num2str(k_eq_delta_I*1e6) '+-' num2str(sigma_k_eq_delta_I*1e6) ' N/m'])

[v1, dv1, sig]=round_significance(k_eq_delta_I*1e6, sigma_k_eq_delta_I*1e6);
fprintf(equipdx_tex,'%s\n',['\newcommand{\kappaequiExpIdelta}{' v1 '\pm' dv1 '}']);



load('Data_positions_Fig9_1P4_S.mat'); 

[k_eq_delta_II,sigma_k_eq_delta_II]=eq1d(x(1:subs:size(x,1),:),T,1e-9);

disp('................')
disp('Equipartition analysis considering deltax=1nm')

disp(['k_eq: ' num2str(k_eq_delta_II*1e6) '+-' num2str(sigma_k_eq_delta_II*1e6) ' N/m'])

[v1, dv1, sig]=round_significance(k_eq_delta_II*1e6, sigma_k_eq_delta_II*1e6);
fprintf(equipdx_tex,'%s\n',['\newcommand{\kappaequiExpIIdelta}{' v1 '\pm' dv1 '}']);


load('Data_positions_Fig9_1P6_S.mat'); 

[k_eq_delta_III,sigma_k_eq_delta_III]=eq1d(x(1:subs:size(x,1),:),T,1e-9);

disp('................')
disp('Equipartition analysis considering deltax=1nm')

disp(['k_eq: ' num2str(k_eq_delta_III*1e6) '+-' num2str(sigma_k_eq_delta_III*1e6) ' N/m'])

[v1, dv1, sig]=round_significance(k_eq_delta_III*1e6, sigma_k_eq_delta_III*1e6);
fprintf(equipdx_tex,'%s\n',['\newcommand{\kappaequiExpIIIdelta}{' v1 '\pm' dv1 '}']);

