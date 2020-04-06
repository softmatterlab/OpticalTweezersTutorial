 
clear all 
close all
%load data files
addpath ../data/
addpath ../statistics_func/
%%
%first experiment
disp('first experiment')
load('Data_x_positions_Exp_I.mat')

 k_th=13*1e-6;



x = x - repmat(mean(x),size(x,1),1);

xl = reshape(x,[size(x,1)*size(x,2),1]);
N=length(xl);

nsubs=3;
subs=1;


[k_bay_s, sigma_k_bay_s, gamma_bay_s, sigma_gamma_bay_s, D_bay_s, sigma_D_bay_s]= bayesian1d(x(:, 1),dt*subs,T, a, nsubs, k_th);

disp('................')

k=k_bay_s*1e6;
dk=sigma_k_bay_s*1e6;
D=D_bay_s*1e12;
dD=sigma_D_bay_s*1e12;
g=1e9*gamma_bay_s;
dg=1e9*sigma_gamma_bay_s;

[v1, dv1, ~]=round_significance(k, dk);


[v1, dv1, ~]=round_significance(D, dD);


[v1, dv1, ~]=round_significance(g, dg);





disp('Bayesian Inference analysis for a single experiment')
 
disp(['k_bay: ' num2str(k_bay_s*1e6) '+-' num2str(sigma_k_bay_s*1e6) ' pN/um'])
 
disp(['D_bay: ' num2str(D_bay_s*1e12) '+-' num2str(sigma_D_bay_s*1e12) ' um^2/s'])
 
disp(['gamma_bay:' num2str(gamma_bay_s*1e9) '+-'  num2str(sigma_gamma_bay_s*1e9) ' pN ms/um ']);
 
disp('................')




Nexp=5;
for jj=1:Nexp
    
[k_bay_ss(jj), sigma_k_bay_ss(jj), gamma_bay_ss(jj), sigma_gamma_bay_ss(jj), D_bay_ss(jj), sigma_D_bay_ss(jj)]= bayesian1d(x(:, jj),dt*subs,T, a, nsubs, k_th);


end

disp('................')

disp('Bayesian Inference analysis averaging all the experiments')
 
disp(['k_bay: ' num2str(mean(k_bay_ss)*1e6) '+-' num2str(std(k_bay_ss)*1e6) ' pN/um'])
 
disp(['D_bay: ' num2str(mean(D_bay_ss)*1e12) '+-' num2str(std(D_bay_ss)*1e12) ' um^2/s'])
 
disp(['gamma_bay:' num2str(mean(gamma_bay_ss)*1e9) '+-'  num2str(std(gamma_bay_ss)*1e9) ' pN ms/um ']);
 
disp('................')


[k_bay_l, sigma_k_bay_l, gamma_bay_l, sigma_gamma_bay_l, D_bay_l, sigma_D_bay_l]= bayesian1d(xl(1:subs:end),dt*subs,T, a, nsubs, k_th);

disp('................')

disp('Bayesian Inference analysis using all the experiments as a big experiment ')
 
disp(['k_bay: ' num2str(k_bay_l*1e6) '+-' num2str(sigma_k_bay_l*1e6) ' pN/um'])
 
disp(['D_bay: ' num2str(D_bay_l*1e12) '+-' num2str(sigma_D_bay_l*1e12) ' um^2/s'])
 
disp(['gamma_bay:' num2str(gamma_bay_l*1e9) '+-'  num2str(sigma_gamma_bay_l*1e9) ' pN ms/um ']);
 
disp('................')


%%
%second experiment
disp('second experiment')


load('Data_x_positions_Exp_II.mat')

k_th=35*1e-6;


x = x - repmat(mean(x),size(x,1),1);

xl = reshape(x,[size(x,1)*size(x,2),1]);
N=length(xl);

nsubs=3;
subs=1;


[k_bay_s, sigma_k_bay_s, gamma_bay_s, sigma_gamma_bay_s, D_bay_s, sigma_D_bay_s]= bayesian1d(x(:, 1),dt*subs,T, a, nsubs, k_th);
disp('................')

k=k_bay_s*1e6;
dk=sigma_k_bay_s*1e6;
D=D_bay_s*1e12;
dD=sigma_D_bay_s*1e12;
g=1e9*gamma_bay_s;
dg=1e9*sigma_gamma_bay_s;

[v1, dv1,~]=round_significance(k, dk);


[v1, dv1,~]=round_significance(D, dD);


[v1, dv1, ~]=round_significance(g, dg);



disp('Bayesian Inference analysis for a single experiment')
 
disp(['k_bay: ' num2str(k_bay_s*1e6) '+-' num2str(sigma_k_bay_s*1e6) ' pN/um'])
 
disp(['D_bay: ' num2str(D_bay_s*1e12) '+-' num2str(sigma_D_bay_s*1e12) ' um^2/s'])
 
disp(['gamma_bay:' num2str(gamma_bay_s*1e9) '+-'  num2str(sigma_gamma_bay_s*1e9) ' pN ms/um ']);
 
disp('................')




Nexp=5;
for jj=1:Nexp
    
[k_bay_ss(jj), sigma_k_bay_ss(jj), gamma_bay_ss(jj), sigma_gamma_bay_ss(jj), D_bay_ss(jj), sigma_D_bay_ss(jj)]= bayesian1d(x(:, jj),dt*subs,T, a, nsubs, k_th);
end

disp('................')

disp('Bayesian Inference analysis averaging all the experiments')
 
disp(['k_bay: ' num2str(mean(k_bay_ss)*1e6) '+-' num2str(std(k_bay_ss)*1e6) ' pN/um'])
 
disp(['D_bay: ' num2str(mean(D_bay_ss)*1e12) '+-' num2str(std(D_bay_ss)*1e12) ' um^2/s'])
 
disp(['gamma_bay:' num2str(mean(gamma_bay_ss)*1e9) '+-'  num2str(std(gamma_bay_ss)*1e9) ' pN ms/um ']);
 
disp('................')


[k_bay_l, sigma_k_bay_l, gamma_bay_l, sigma_gamma_bay_l, D_bay_l, sigma_D_bay_l]= bayesian1d(xl(1:subs:end),dt*subs,T, a, nsubs, k_th);
disp('................')

disp('Bayesian Inference analysis using all the experiments as a big experiment ')
 
disp(['k_bay: ' num2str(k_bay_l*1e6) '+-' num2str(sigma_k_bay_l*1e6) ' pN/um'])
 
disp(['D_bay: ' num2str(D_bay_l*1e12) '+-' num2str(sigma_D_bay_l*1e12) ' um^2/s'])
 
disp(['gamma_bay:' num2str(gamma_bay_l*1e9) '+-'  num2str(sigma_gamma_bay_l*1e9) ' pN ms/um ']);
 
disp('................')


%%
%third experiment
disp('third experiment')
load('Data_x_positions_Exp_III.mat')

k_th=58*1e-6;

x = x - repmat(mean(x),size(x,1),1);

xl = reshape(x,[size(x,1)*size(x,2),1]);
N=length(xl);

nsubs=3;
subs=1;


[k_bay_s, sigma_k_bay_s, gamma_bay_s, sigma_gamma_bay_s, D_bay_s, sigma_D_bay_s]= bayesian1d(x(:, 1),dt*subs,T, a, nsubs, k_th);
disp('................')

k=k_bay_s*1e6;
dk=sigma_k_bay_s*1e6;
D=D_bay_s*1e12;
dD=sigma_D_bay_s*1e12;
g=1e9*gamma_bay_s;
dg=1e9*sigma_gamma_bay_s;

[v1, dv1,~]=round_significance(k, dk);


[v1, dv1, ~]=round_significance(D, dD);


[v1, dv1, ~]=round_significance(g, dg);



disp('Bayesian Inference analysis for a single experiment')
 
disp(['k_bay: ' num2str(k_bay_s*1e6) '+-' num2str(sigma_k_bay_s*1e6) ' pN/um'])
 
disp(['D_bay: ' num2str(D_bay_s*1e12) '+-' num2str(sigma_D_bay_s*1e12) ' um^2/s'])
 
disp(['gamma_bay:' num2str(gamma_bay_s*1e9) '+-'  num2str(sigma_gamma_bay_s*1e9) ' pN ms/um ']);
 
disp('................')




Nexp=5;
for jj=1:Nexp
    
[k_bay_ss(jj), sigma_k_bay_ss(jj), gamma_bay_ss(jj), sigma_gamma_bay_ss(jj), D_bay_ss(jj), sigma_D_bay_ss(jj)]= bayesian1d(x(:, jj),dt*subs,T, a, nsubs, k_th);

end

disp('................')

disp('Bayesian Inference analysis averaging all the experiments')
 
disp(['k_bay: ' num2str(mean(k_bay_ss)*1e6) '+-' num2str(std(k_bay_ss)*1e6) ' pN/um'])
 
disp(['D_bay: ' num2str(mean(D_bay_ss)*1e12) '+-' num2str(std(D_bay_ss)*1e12) ' um^2/s'])
 
disp(['gamma_bay:' num2str(mean(gamma_bay_ss)*1e9) '+-'  num2str(std(gamma_bay_ss)*1e9) ' pN ms/um ']);
 
disp('................')


[k_bay_l, sigma_k_bay_l, gamma_bay_l, sigma_gamma_bay_l, D_bay_l, sigma_D_bay_l]= bayesian1d(xl(1:subs:end),dt*subs,T, a, nsubs, k_th);
disp('................')

disp('Bayesian Inference analysis using all the experiments as a big experiment ')
 
disp(['k_bay: ' num2str(k_bay_l*1e6) '+-' num2str(sigma_k_bay_l*1e6) ' pN/um'])
 
disp(['D_bay: ' num2str(D_bay_l*1e12) '+-' num2str(sigma_D_bay_l*1e12) ' um^2/s'])
 
disp(['gamma_bay:' num2str(gamma_bay_l*1e9) '+-'  num2str(sigma_gamma_bay_l*1e9) ' pN ms/um ']);
 
disp('................')



