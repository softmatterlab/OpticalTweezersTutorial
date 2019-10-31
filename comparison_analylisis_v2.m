clear all
close all
addpath pot
addpath eq
addpath acf
addpath psd
addpath msd
addpath forma
addpath bayesian
load('Data_positions_Fig9_1P2_S.mat')



x = x - repmat(mean(x),size(x,1),1);
kB=1.38e-23; % Boltzmann constant [m^2kg/s^2K]
T=300;  % Temperature [K]
r=1.03E-6;      % Particle radius [m]
v=0.00002414*10^(247.8/(-140+T));  % Water viscosity [Pa*s]
gamma=pi*6*r*v; %[m*Pa*s]
dnn=5e3:2.5*1e2:1e5;
disp(length(dnn));
subs=10;
xs=x(1:subs:end, :);

%%
for nj=1:length(dnn)


disp(nj);
disp(dnn(nj));
xn =xs(1:dnn(nj), :);
N=size(xn);



%Potential
[k_pot(nj), sigma_k_pot(nj), ~, ~, ~, ~, ~, ~, ~,~]=pot_lfit(xn,T,50);
%Equipartition
[k_eq(nj), sigma_k_eq(nj)]=eq1d(xn,T);
%autocorrelation function
[k_acf(nj), sigma_k_acf(nj), D_acf(nj), sigma_D_acf(nj),gamma_exp(nj), sigma_gamma_exp(nj),~, ~, ~,~, ~, ~]=acf_lfit(xn,T,dt*subs);
end

%%

%power spectral density
subs=10;
nw=round(size(xs,1)/500); 
for nj=1:length(dnn)
disp(nj);
disp(dnn(nj));
xn =xs(1:dnn(nj), :);
N=size(xn);
% number of windows

[mfc_psd,D_psd(nj),Efc_psd,sigma_D_psd,f,XX,fw_mean,Pk,EPk,fcut]=psd_lfit(xn,dt*subs,nw,1/4);

mgamma_psd=kB*T./D_psd(nj);

sigma_gamma_psd=kB*T./D_psd(nj)^2*sigma_D_psd;

k2_psd(nj)=2*pi*gamma*mfc_psd;

sigma_k_psd=2*pi*gamma*Efc_psd;

% estimation of k using the estimated gamma
k_psd(nj)=2*pi*mgamma_psd.*mfc_psd;

sigma_k2_psd(nj)=2*pi*mgamma_psd.*Efc_psd+2*pi*mfc_psd*sigma_gamma_psd;
end



   
%%
%especial data format for FORMA and bayesian  analysis
subs=10;
xr=reshape(x, [size(x,1)*size(x,2),1 ]);
xr=xr(1:subs:end)
for nj=1:length(dnn)
disp(nj);
disp(dnn(nj));
xn =xr(1:dnn(nj), :);
N=size(xn);

%FORMA

[fc_forma,D_forv,Efc_forma,sigma_D_forma(nj)] = forma1d(xn, dt*subs);
D_forma(nj)=mean(D_forv);
gamma_forma=kB*T./D_forma(nj);
kfor=gamma_forma.*mean(fc_forma);

k_forma(nj)=kfor;
sigma_k_forma(nj)=gamma_forma.*Efc_forma;
%BAYESIAN
[k_bay(nj), sigma_k_bay(nj), gamma_bay(nj), sigma_gamma_bay(nj), D_bay(nj), sigma_D_bay(nj)]= bayesian(xn, dt*subs,T, a);

end


%%
for nj=1:length(dnn)
 xn =xs(1:dnn(nj), :);
%mean square displacement
maxlag=50;
[k_msd(nj),sigma_k_msd(nj),~, ~, D_msd(nj), sigma_D_msd(nj), ~, ~, ~, ~, gamma_msd(nj), sigma2_gamma_msd(nj)]=msd_nfilt(xn,T,dt*subs,maxlag);


end
%%
figure(1)
%load('results_comparison.mat')
semilogx(dnn, k_pot*1e6, 'DisplayName', 'Potential', 'Color', 'red', 'LineWidth', 3)
hold on
semilogx(dnn, k_eq*1e6, 'DisplayName', 'Equipartition', 'Color', 'blue', 'LineWidth', 3)
semilogx(dnn, k_acf*1e6, 'DisplayName', 'ACF', 'Color', 'cyan', 'LineWidth', 3)
semilogx(dnn, k_psd*1e6, 'DisplayName', 'PSD', 'Color', 'yellow', 'LineWidth', 3)
semilogx(dnn, k_msd*1e6, 'DisplayName', 'MSD', 'Color', 'black', 'LineWidth', 3)

semilogx(dnn, k_forma*1e6, 'DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 3)
hold on
semilogx(dnn, k_bay*1e6, 'DisplayName', 'BAYESIAN', 'Color', [0.1 0.5 0.7], 'LineWidth', 3)
legend
xlabel('$N_s$', 'Interpreter', 'Latex', 'FontSize',30)
ylabel('$k( p \rm N  \mu \rm m^{-1})$', 'Interpreter', 'Latex', 'FontSize',30)
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
hold off
%%
%load('results_comparison.mat')
figure(2)
semilogx(dnn, sigma_k_pot./k_pot, 'DisplayName', 'Potential', 'Color', 'red', 'LineWidth', 3)
hold on
semilogx(dnn, sigma_k_eq./k_eq, 'DisplayName', 'Equipartition', 'Color', 'blue', 'LineWidth', 3)
semilogx(dnn, sigma_k_acf./k_acf, 'DisplayName', 'ACF', 'Color', 'cyan', 'LineWidth', 3)
semilogx(dnn, sigma_k_psd./k_psd, 'DisplayName', 'PSD', 'Color', 'yellow', 'LineWidth', 3)
semilogx(dnn, sigma_k_msd./k_msd, 'DisplayName', 'MSD', 'Color', 'black', 'LineWidth', 3)

semilogx(dnn, sigma_k_forma./k_forma, 'DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 3)

semilogx(dnn, sigma_k_bay./k_bay, 'DisplayName', 'BAYESIAN', 'Color', [0.1 0.5 0.7], 'LineWidth', 3)
legend
xlabel('$N_s$', 'Interpreter', 'Latex', 'FontSize',30)
ylabel('$\sigma k/k$', 'Interpreter', 'Latex', 'FontSize',30)
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
hold off

%%
figure(3)


semilogx(dnn, D_acf, 'DisplayName', 'ACF', 'Color', 'cyan', 'LineWidth', 3)
hold on
semilogx(dnn, D_psd*1e6, 'DisplayName', 'PSD', 'Color', 'yellow', 'LineWidth', 3)
semilogx(dnn, D_msd*1e6, 'DisplayName', 'MSD', 'Color', 'black', 'LineWidth', 3)

semilogx(dnn, D_forma*1e6, 'DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 3)

semilogx(dnn, D_bay*1e6, 'DisplayName', 'BAYESIAN', 'Color', [0.1 0.5 0.7], 'LineWidth', 3)
legend
hold off
%%

%load('results_comparison.mat')
figure(4)

%semilogx(dnn, D_eq*1e6, 'DisplayName', 'Equipartition', 'Color', 'blue', 'LineWidth', 3)
semilogx(dnn, D_acf./sigma_D_acf, 'DisplayName', 'ACF', 'Color', 'cyan', 'LineWidth', 3)
hold on
semilogx(dnn, D_psd./sigma_D_psd, 'DisplayName', 'PSD', 'Color', 'yellow', 'LineWidth', 3)
semilogx(dnn, D_msd./sigma_D_msd, 'DisplayName', 'MSD', 'Color', 'black', 'LineWidth', 3)
semilogx(dnn, D_forma./sigma_D_forma, 'DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 3)

semilogx(dnn, D_bay./sigma_D_bay, 'DisplayName', 'BAYESIAN', 'Color', [0.1 0.5 0.7], 'LineWidth', 3)
legend
hold off