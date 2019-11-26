clear all
close all
addpath pot
addpath eq
addpath acf
addpath psd
addpath msd
addpath forma
addpath bayesian
filename='Data_positions_Fig9_1P2_S.mat';
load(filename)



x = x - repmat(mean(x),size(x,1),1);
kB=1.38064852e-23; % Boltzmann constant [m^2kg/s^2K]

[Nfull, Nexp]=size(x);
v=0.00002414*10^(247.8/(-140+T));  % Water viscosity [Pa*s]
gamma=pi*6*a*v; %[m*Pa*s]

subs=1;
dnn=floor(1e1:5*5e2:1e6);
disp(size(dnn));


%%
for nns=1:length(dnn)-1


for nexp=1:Nexp
    %first we create the Nexp different experiments
xs=x(1:dnn(nns),nexp);

nsubs=3;

%FORMA

[fc_forma,D_forv,Efc_forma(nexp, nns+1),sigma_D_forma(nexp, nns+1)] = forma1d(xs, dt, nsubs);
D_forma(nexp, nns+1)=mean(D_forv);
gamma_forma(nexp, nns+1)=kB*T./D_forma(nexp, nns+1);
kfor=gamma_forma(nexp, nns+1).*mean(fc_forma);

k_forma(nexp, nns+1)=kfor;
sigma_k_forma(nexp, nns+1)=gamma_forma(nexp, nns+1).*Efc_forma(nexp, nns+1);
%BAYESIAN
[k_bay(nexp, nns+1), sigma_k_bay(nexp, nns+1), gamma_bay(nexp, nns+1), sigma_gamma_bay(nexp, nns+1), D_bay(nexp, nns+1), sigma_D_bay(nexp, nns+1)]= bayesian(xs, dt,T, a, nsubs , 14*1e-6);
end


end

%%


kk_forma=mean(k_forma, 1);
sigma_kk_forma=std(k_forma,[], 1);

DD_forma=mean(D_forma, 1);
sigma_DD_forma=std(D_forma,[], 1);

kk_bay=mean(k_bay,1);
sigma_kk_bay=std(k_bay, [] , 1);

DD_bay=mean(D_bay, 1);
sigma_DD_bay=std(D_bay,[], 1);
%Average over experiments for each number of points in the trajectory


%%




% save(['results_averaged_bayes_vs_forma' , filename ,datestr(now, 'yyyymmddTHHMMSS'), '.mat'],'subs', 'dnn',...
%     'kk_forma', 'sigma_kk_forma', 'kk_bay', 'sigma_kk_bay','DD_forma',...
%     'sigma_DD_forma', 'DD_bay', 'sigma_DD_bay', 'gg_forma', 'sigma_gg_forma', 'gg_bay', 'sigma_gg_bay')
% %%
%dnn=dnn(1:50);
figure(1)
%load('results_comparison.mat')
% %semilogx(dnn, kk_pot, 'DisplayName', 'Potential', 'Color', 'red', 'LineWidth', 3)
% %hold on
% semilogx(dnn, kk_eq, 'DisplayName', 'Equipartition', 'Color', 'blue', 'LineWidth', 3)
% semilogx(dnn, kk_acf, 'DisplayName', 'ACF', 'Color', 'cyan', 'LineWidth', 3)
% semilogx(dnn, kk_psd, 'DisplayName', 'PSD', 'Color', 'yellow', 'LineWidth', 3)
% semilogx(dnn, kk_msd, 'DisplayName', 'MSD', 'Color', 'black', 'LineWidth', 3)

semilogx(dnn, kk_forma, 'DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 3)
hold on
semilogx(dnn, kk_bay, '--','DisplayName', 'BAYESIAN', 'Color', [0.1 0.5 0.7], 'LineWidth', 3)
legend
xlabel('$N_s$', 'Interpreter', 'Latex', 'FontSize',30)
ylabel('$k( p \rm N  \mu \rm m^{-1})$', 'Interpreter', 'Latex', 'FontSize',30)
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
hold off
%%
%load('results_comparison.mat')
figure(2)
%semilogx(dnn, sigma_kk_pot./kk_pot, 'DisplayName', 'Potential', 'Color', 'red', 'LineWidth', 3)

%semilogx(dnn, sigma_kk_eq./kk_eq, 'DisplayName', 'Equipartition', 'Color', 'blue', 'LineWidth', 3)
%semilogx(dnn, sigma_kk_acf./kk_acf, 'DisplayName', 'ACF', 'Color', 'cyan', 'LineWidth', 3)
%semilogx(dnn, sigma_kk_psd./kk_psd, 'DisplayName', 'PSD', 'Color', 'yellow', 'LineWidth', 3)
%semilogx(dnn, sigma_kk_msd./kk_msd, 'DisplayName', 'MSD', 'Color', 'black', 'LineWidth', 3)


semilogx(dnn, sigma_kk_forma./kk_forma, 'DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 3)

hold on
semilogx(dnn, sigma_kk_bay./kk_bay , '--','DisplayName', 'BAYESIAN',  'Color', [0.1 0.5 0.7], 'LineWidth', 3)
legend
xlabel('$N_s$', 'Interpreter', 'Latex', 'FontSize',30)
ylabel('$\Delta k/k$', 'Interpreter', 'Latex', 'FontSize',30)
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
hold off

%%
figure(3)


% semilogx(dnn, DD_acf*1e6, 'DisplayName', 'ACF', 'Color', 'cyan', 'LineWidth', 3)
% hold on
% semilogx(dnn, DD_psd*1e6, 'DisplayName', 'PSD', 'Color', 'yellow', 'LineWidth', 3)
% semilogx(dnn, DD_msd*1e6, 'DisplayName', 'MSD', 'Color', 'black', 'LineWidth', 3)

semilogx(dnn, DD_forma*1e6, 'DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 3)
hold on
semilogx(dnn, DD_bay*1e6,'--','DisplayName', 'BAYESIAN', 'Color', [0.1 0.5 0.7], 'LineWidth', 3)
legend
xlabel('$N_s$', 'Interpreter', 'Latex', 'FontSize',30)
ylabel('$D(\mu m^2/s)$', 'Interpreter', 'Latex', 'FontSize',30)
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
hold off

%%

%load('results_comparison.mat')
figure(4)

%semilogx(dnn, D_eq*1e6, 'DisplayName', 'Equipartition', 'Color', 'blue', 'LineWidth', 3)
% semilogx(dnn,sigma_DD_acf./ DD_acf, 'DisplayName', 'ACF', 'Color', 'cyan', 'LineWidth', 3)
% hold on
% semilogx(dnn, sigma_DD_psd./DD_psd, 'DisplayName', 'PSD', 'Color', 'yellow', 'LineWidth', 3)
% semilogx(dnn, sigma_DD_msd./DD_msd, 'DisplayName', 'MSD', 'Color', 'black', 'LineWidth', 3)
semilogx(dnn, sigma_DD_forma./DD_forma, 'DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 3)
hold on
semilogx(dnn,sigma_DD_bay./DD_bay,'--', 'DisplayName', 'BAYESIAN', 'Color', [0.1 0.5 0.7], 'LineWidth', 3)
legend
xlabel('$N_s$', 'Interpreter', 'Latex', 'FontSize',30)
ylabel('$\Delta D/D$', 'Interpreter', 'Latex', 'FontSize',30)
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
hold off


% figure(5)
% %load('results_comparison.mat')
% 
% % semilogx(dnn, gg_acf, 'DisplayName', 'ACF', 'Color', 'cyan', 'LineWidth', 3)
% % hold on
% % semilogx(dnn, gg_psd, 'DisplayName', 'PSD', 'Color', 'yellow', 'LineWidth', 3)
% % semilogx(dnn, gg_msd, 'DisplayName', 'MSD', 'Color', 'black', 'LineWidth', 3)
% 
% semilogx(dnn, gg_forma, 'DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 3)
% hold on
% semilogx(dnn, gg_bay, 'DisplayName', 'BAYESIAN', 'Color', [0.1 0.5 0.7], 'LineWidth', 3)
% legend
% xlabel('$N_s$', 'Interpreter', 'Latex', 'FontSize',30)
% ylabel('$k( p \rm N  \mu \rm m^{-1})$', 'Interpreter', 'Latex', 'FontSize',30)
% set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
% hold off