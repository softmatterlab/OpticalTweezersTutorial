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
kB=1.38e-23; % Boltzmann constant [m^2kg/s^2K]
T=300;  % Temperature [K]
r=1.03E-6;      % Particle radius [m]
v=0.00002414*10^(247.8/(-140+T));  % Water viscosity [Pa*s]
gamma=pi*6*r*v; %[m*Pa*s]
dnn=5e3:2.5*1e2:1e5;
disp(length(dnn));
subs=10;
for nexp=1:10
xs=x(nexp:subs:end, :);

%%
for nj=1:length(dnn)


disp(nj);
disp(dnn(nj));

xn =xs(1:dnn(nj), :);
N=size(xn);



%Potential
[k_pot(nj, nexp), sigma_k_pot(nj, nexp), ~, ~, ~, ~, ~, ~, ~,~]=pot_lfit(xn,T,100);
%Equipartition
[k_eq(nj, nexp), sigma_k_eq(nj, nexp)]=eq1d(xn,T);
%autocorrelation function
[k_acf(nj, nexp), sigma_k_acf(nj, nexp), D_acf(nj, nexp), sigma_D_acf(nj, nexp),gamma_exp(nj, nexp), sigma_gamma_exp(nj, nexp),~, ~, ~,~, ~, ~]=acf_lfit(xn,T,dt*subs);
end

%%

%power spectral density

nw=round(size(xs,1)/500); 
for nj=1:length(dnn)
disp(nj);
disp(dnn(nj));
xn =xs(1:dnn(nj), :);
N=size(xn);
% number of windows

[mfc_psd,D_psd(nj, nexp),Efc_psd,sigma_D_psd,f,XX,fw_mean,Pk,EPk,fcut]=psd_lfit(xn,dt*subs,nw,1/4);

mgamma_psd=kB*T./D_psd(nj, nexp);

sigma_gamma_psd=kB*T./D_psd(nj)^2*sigma_D_psd;

k2_psd(nj, nexp)=2*pi*gamma*mfc_psd;

sigma_k_psd=2*pi*gamma*Efc_psd;

% estimation of k using the estimated gamma
k_psd(nj, nexp)=2*pi*mgamma_psd.*mfc_psd;

sigma_k2_psd(nj, nexp)=2*pi*mgamma_psd.*Efc_psd+2*pi*mfc_psd*sigma_gamma_psd;
end



   
%%
%especial data format for FORMA and bayesian  analysis

xr=reshape(xs, [size(xs,1)*size(xs,2),1 ]);
%xr=xr(1:subs:end)
for nj=1:length(dnn)
disp(nj);
disp(dnn(nj));
xn =xr(1:dnn(nj), :);
N=size(xn);

%FORMA

[fc_forma,D_forv,Efc_forma,sigma_D_forma(nj, nexp)] = forma1d(xn, dt*subs);
D_forma(nj, nexp)=mean(D_forv);
gamma_forma=kB*T./D_forma(nj);
kfor=gamma_forma.*mean(fc_forma);

k_forma(nj, nexp)=kfor;
sigma_k_forma(nj, nexp)=gamma_forma.*Efc_forma;
%BAYESIAN
[k_bay(nj, nexp), sigma_k_bay(nj, nexp), gamma_bay(nj, nexp), sigma_gamma_bay(nj, nexp), D_bay(nj, nexp), sigma_D_bay(nj, nexp)]= bayesian(xn, dt*subs,T, a);

end


%%
for nj=1:length(dnn)
 xn =xs(1:dnn(nj), :);
%mean square displacement
maxlag=50;
[k_msd(nj, nexp),sigma_k_msd(nj, nexp),~, ~, D_msd(nj, nexp), sigma_D_msd(nj, nexp), ~, ~, ~, ~, gamma_msd(nj, nexp), sigma2_gamma_msd(nj, nexp)]=msd_nfilt(xn,T,dt*subs,maxlag);


end

end
%%



for nj=1:length(dnn)
    
kk_pot(nj)=mean(k_pot(nj, :));
sigma_kk_pot(nj)=std(k_pot(nj, :));

kk_eq(nj)=mean(k_eq(nj, :));
sigma_kk_eq(nj)=std(k_eq(nj, :));


kk_acf(nj)=mean(k_acf(nj, :));
sigma_kk_acf(nj)=std(k_acf(nj, :));


kk_psd(nj)=mean(k_psd(nj, :));
sigma_kk_psd(nj)=std(k_psd(nj, :));


kk_forma(nj)=mean(k_forma(nj, :));
sigma_kk_forma(nj)=std(k_forma(nj, :));


kk_bay(nj)=mean(k_bay(nj, :));
sigma_kk_bay(nj)=std(k_bay(nj, :));


kk_msd(nj)=mean(k_msd(nj, :));
sigma_kk_msd(nj)=std(k_msd(nj, :));
%%Diffusion coeficcient

DD_msd(nj)=mean(D_msd(nj, :));
sigma_DD_msd(nj)=std(D_msd(nj, :));



DD_acf(nj)=mean(D_acf(nj, :));
sigma_DD_acf(nj)=std(D_acf(nj, :));


DD_psd(nj)=mean(D_psd(nj, :));
sigma_DD_psd(nj)=std(D_psd(nj, :));


DD_forma(nj)=mean(D_forma(nj, :));
sigma_DD_forma(nj)=std(D_forma(nj, :));


DD_bay(nj)=mean(D_bay(nj, :));
sigma_DD_bay(nj)=std(D_bay(nj, :));

DD_msd(nj)=mean(D_msd(nj, :));
sigma_DD_msd(nj)=std(D_msd(nj, :));


end
save(['results_averaged_' , filename ,datestr(now, 'yyyymmddTHHMMSS'), '.mat'],'subs', 'dnn', 'kk_pot', 'sigma_kk_pot',...
    'kk_eq', 'sigma_kk_eq','kk_acf', 'sigma_kk_acf', 'kk_psd', 'sigma_kk_psd',...
    'kk_msd', 'sigma_kk_msd', 'kk_forma', 'sigma_kk_forma', 'kk_bay', 'sigma_kk_bay', ...
    'DD_acf', 'sigma_DD_acf', 'DD_psd', 'sigma_DD_psd', 'DD_msd', 'sigma_DD_msd','DD_forma',...
    'sigma_DD_forma', 'DD_bay', 'sigma_DD_bay')
%%


figure(1)
%load('results_comparison.mat')
semilogx(dnn, kk_pot, 'DisplayName', 'Potential', 'Color', 'red', 'LineWidth', 3)
hold on
semilogx(dnn, kk_eq, 'DisplayName', 'Equipartition', 'Color', 'blue', 'LineWidth', 3)
semilogx(dnn, kk_acf, 'DisplayName', 'ACF', 'Color', 'cyan', 'LineWidth', 3)
semilogx(dnn, kk_psd, 'DisplayName', 'PSD', 'Color', 'yellow', 'LineWidth', 3)
semilogx(dnn, kk_msd, 'DisplayName', 'MSD', 'Color', 'black', 'LineWidth', 3)

semilogx(dnn, kk_forma, 'DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 3)
hold on
semilogx(dnn, kk_bay, 'DisplayName', 'BAYESIAN', 'Color', [0.1 0.5 0.7], 'LineWidth', 3)
legend
xlabel('$N_s$', 'Interpreter', 'Latex', 'FontSize',30)
ylabel('$k( p \rm N  \mu \rm m^{-1})$', 'Interpreter', 'Latex', 'FontSize',30)
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
hold off
%%
%load('results_comparison.mat')
figure(2)
semilogx(dnn, sigma_kk_pot./kk_pot, 'DisplayName', 'Potential', 'Color', 'red', 'LineWidth', 3)
hold on
semilogx(dnn, sigma_kk_eq./kk_eq, 'DisplayName', 'Equipartition', 'Color', 'blue', 'LineWidth', 3)
semilogx(dnn, sigma_kk_acf./kk_acf, 'DisplayName', 'ACF', 'Color', 'cyan', 'LineWidth', 3)
semilogx(dnn, sigma_kk_psd./kk_psd, 'DisplayName', 'PSD', 'Color', 'yellow', 'LineWidth', 3)
semilogx(dnn, sigma_kk_msd./kk_msd, 'DisplayName', 'MSD', 'Color', 'black', 'LineWidth', 3)


semilogx(dnn, sigma_kk_forma./kk_forma, 'DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 3)

semilogx(dnn, sigma_kk_bay./kk_bay ,'DisplayName', 'BAYESIAN',  'Color', [0.1 0.5 0.7], 'LineWidth', 3)
legend
xlabel('$N_s$', 'Interpreter', 'Latex', 'FontSize',30)
ylabel('$\Delta k/k$', 'Interpreter', 'Latex', 'FontSize',30)
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
hold off

%%
figure(3)


semilogx(dnn, DD_acf*1e6, 'DisplayName', 'ACF', 'Color', 'cyan', 'LineWidth', 3)
hold on
semilogx(dnn, DD_psd*1e6, 'DisplayName', 'PSD', 'Color', 'yellow', 'LineWidth', 3)
semilogx(dnn, DD_msd*1e6, 'DisplayName', 'MSD', 'Color', 'black', 'LineWidth', 3)

semilogx(dnn, DD_forma*1e6, 'DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 3)

semilogx(dnn, DD_bay*1e6, 'DisplayName', 'BAYESIAN', 'Color', [0.1 0.5 0.7], 'LineWidth', 3)
legend
xlabel('$N_s$', 'Interpreter', 'Latex', 'FontSize',30)
ylabel('$D(\mu m^2/s)$', 'Interpreter', 'Latex', 'FontSize',30)
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
hold off

%%

%load('results_comparison.mat')
figure(4)

%semilogx(dnn, D_eq*1e6, 'DisplayName', 'Equipartition', 'Color', 'blue', 'LineWidth', 3)
semilogx(dnn,sigma_DD_acf./ DD_acf, 'DisplayName', 'ACF', 'Color', 'cyan', 'LineWidth', 3)
hold on
semilogx(dnn, sigma_DD_psd./DD_psd, 'DisplayName', 'PSD', 'Color', 'yellow', 'LineWidth', 3)
semilogx(dnn, sigma_DD_msd./DD_msd, 'DisplayName', 'MSD', 'Color', 'black', 'LineWidth', 3)
semilogx(dnn, sigma_DD_forma./DD_forma, 'DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 3)

semilogx(dnn,sigma_DD_bay./DD_bay, 'DisplayName', 'BAYESIAN', 'Color', [0.1 0.5 0.7], 'LineWidth', 3)
legend
xlabel('$N_s$', 'Interpreter', 'Latex', 'FontSize',30)
ylabel('$\Delta D/D$', 'Interpreter', 'Latex', 'FontSize',30)
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
hold off