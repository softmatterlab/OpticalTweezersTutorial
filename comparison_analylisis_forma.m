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




kB=1.38e-23; % Boltzmann constant [m^2kg/s^2K]
T=300;  % Temperature [K]
r=1.03E-6;      % Particle radius [m]
v=0.00002414*10^(247.8/(-140+T));  % Water viscosity [Pa*s]
gamma=pi*6*r*v; %[m*Pa*s]
dnn=5e2*(1:1e2:2000-1);
disp(length(dnn));

for nj=1:length(dnn)

    
disp(nj);
xn =x(1:dnn(nj)-1, :);

dxn=diff(x(1:dnn(nj), :), 1, 1);

subs=10;

[fc_forma,D_forv,Efc_forma,sigma_D_forma_v] = forma1d_v2(xn(1:subs:end, :), dxn(1:subs:end, :),dt*subs);
D_forma(nj)=mean(D_forv);
gamma_forma=kB*T./D_forma(nj);
k_forma=gamma*fc_forma;
kfor=gamma_forma.*mean(fc_forma);
disp(kfor);
k_forma(nj)=kfor;
sigma_k_forma=gamma_forma.*Efc_forma;



mD_forma=mean(D_forma);

ED_forma=std(D_forma);

gamma_forma=kb*T./D_forma;

mgamma_forma=mean(gamma_forma);

Egamma_forma=std(gamma_forma);

k_forma=gamma*fc_forma;

mk_forma=mean(k_forma);

Ek_forma=std(k_forma);

% estimation of k using the estimated gamma
k2_forma=gamma_forma.*fc_forma;

mk2_forma=mean(k2_forma);

Ek2_forma=std(k2_forma);

Eksol=2*pi*gamma_forma.*Efcsol;
%BAYESIAN
%[k_bay(nj), sigma_k_bay(nj), gamma_bay(nj), sigma_gamma_bay(nj), D_bay(nj), sigma_D_bay(nj)]= bayesian(xl(1:subs:end), dt*subs,T, a);

end


%save('results_comparison.mat', 'k_pot', 'sigma_k_pot', 'k_eq', 'sigma_k_eq')
%%
close all
%load('results_comparison.mat')
% semilogx(dnn(1:30), k_pot*1e6, 'DisplayName', 'Potential', 'Color', 'red', 'LineWidth', 3)
% hold on
% semilogx(dnn(1:30), k_eq*1e6, 'DisplayName', 'Equipartition', 'Color', 'blue', 'LineWidth', 3)
% semilogx(dnn(1:30), k_acf*1e6, 'DisplayName', 'ACF', 'Color', 'cyan', 'LineWidth', 3)
% semilogx(dnn(1:30), k_psd*1e6, 'DisplayName', 'PSD', 'Color', 'yellow', 'LineWidth', 3)
% semilogx(dnn(1:30), k_msd*1e6, 'DisplayName', 'MSD', 'Color', 'black', 'LineWidth', 3)

semilogx(dnn, k_forma, 'DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 3)
hold on
semilogx(dnn, k_bay, 'DisplayName', 'BAYESIAN', 'Color', [0.1 0.5 0.7], 'LineWidth', 3)
legend