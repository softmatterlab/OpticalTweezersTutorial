% Equipartition analysis for Nexp

% Initialization of the workspace
clear;

close all;
addpath ../data/

load('Data_x_positions_Exp_III.mat');
deltax=0;
subs=1; %use a subsampled data set
[k_eq,sigma_k_eq]=eq1d(x(1:subs:size(x,1),:),T,deltax);

disp('................')
disp('Equipartition analysis considering deltax equal to zero')

disp(['k_eq: ' num2str(k_eq*1e6) '+-' num2str(sigma_k_eq*1e6) ' N/m'])

%


load('Data_x_positions_Exp_I.mat'); 

[k_eq_delta_I,sigma_k_eq_delta_I]=eq1d(x(1:subs:size(x,1),:),T,1e-9);

disp('................')
disp('Equipartition analysis considering deltax=1nm')

disp(['k_eq: ' num2str(k_eq_delta_I*1e6) '+-' num2str(sigma_k_eq_delta_I*1e6) ' N/m'])

[v1, dv1, sig]=round_significance(k_eq_delta_I*1e6, sigma_k_eq_delta_I*1e6);


load('Data_x_positions_Exp_II.mat'); 

[k_eq_delta_II,sigma_k_eq_delta_II]=eq1d(x(1:subs:size(x,1),:),T,1e-9);

disp('................')
disp('Equipartition analysis considering deltax=1nm')

disp(['k_eq: ' num2str(k_eq_delta_II*1e6) '+-' num2str(sigma_k_eq_delta_II*1e6) ' N/m'])

[v1, dv1, sig]=round_significance(k_eq_delta_II*1e6, sigma_k_eq_delta_II*1e6);


load('Data_x_positions_Exp_III.mat'); 

[k_eq_delta_III,sigma_k_eq_delta_III]=eq1d(x(1:subs:size(x,1),:),T,1e-9);

disp('................')
disp('Equipartition analysis considering deltax=1nm')

disp(['k_eq: ' num2str(k_eq_delta_III*1e6) '+-' num2str(sigma_k_eq_delta_III*1e6) ' N/m'])

[v1, dv1, sig]=round_significance(k_eq_delta_III*1e6, sigma_k_eq_delta_III*1e6);

