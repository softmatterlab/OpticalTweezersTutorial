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

T=293.15;  % Temperature [K]
r=1.03E-6;      % Particle radius [m]
v=0.00002414*10^(247.8/(-140+T));  % Water viscosity [Pa*s]
gamma=pi*6*r*v; %[m*Pa*s]
subs=3;
dnn=round((1e2:2.5*1e3:1e6)/subs);
disp(length(dnn));

for nexp=1:3
    %first we create the 10 different experiments
xs=x(nexp:subs:end, :);


for Ns=1:length(dnn)
%then we start to change the number of

disp(Ns);
disp(dnn(Ns));

xn =xs(1:dnn(Ns), :);
N=size(xn);
    
%Potential
nbins_pot=50;
try

    [k_pot(Ns, nexp), sigma_k_pot(Ns, nexp), ~, ~, ~, ~, ~, ~, ~,~]=pot_lfit(xn,T,nbins_pot);
catch
    k_pot(Ns, nexp)=NaN;
    sigma_k_pot(Ns, nexp) =NaN;
end

%Equipartition
try
    [k_eq(Ns, nexp), sigma_k_eq(Ns, nexp)]=eq1d(xn,T);
catch
    k_eq(Ns, nexp)=NaN;
    sigma_k_eq(Ns, nexp) =NaN;
end
%autocorrelation function
try
    
    [k_acf(Ns, nexp), sigma_k_acf(Ns, nexp), D_acf(Ns, nexp), sigma_D_acf(Ns, nexp),...
    gamma_exp(Ns, nexp), sigma_gamma_exp(Ns, nexp),~, ~, ~,~, ~, ~]=acf_lfit(xn,T,dt*subs);
catch
    k_acf(Ns, nexp)=NaN;
    sigma_k_acf(Ns, nexp) =NaN;
    D_acf(Ns, nexp)=NaN;
    sigma_D_acf(Ns, nexp)=NaN;
end
%power spectral density
% number of windows
nw_psd=500;
nw=round(1*size(xn,1)/nw_psd); 
NFf=1/4;
try
    [mfc_psd,D_psd(Ns, nexp),Efc_psd,sigma_D_psd,f,XX,fw_mean,Pk,EPk,fcut]=psd_lfit(xn,dt*subs,nw,NFf);
    
    mgamma_psd=kB*T./D_psd(Ns, nexp);
    
    sigma_gamma_psd=kB*T./D_psd(Ns)^2*sigma_D_psd;
    
    k_psd(Ns, nexp)=2*pi*gamma*mfc_psd;
    
    sigma_k_psd=2*pi*gamma*Efc_psd;
    
    % estimation of k using the estimated gamma
    k_psd_var(Ns, nexp)=2*pi*mgamma_psd.*mfc_psd;
    
    sigma_k_psd_var(Ns, nexp)=2*pi*mgamma_psd.*Efc_psd+2*pi*mfc_psd*sigma_gamma_psd;
catch
    k_psd(Ns, nexp)=NaN;
    sigma_k_psd(Ns, nexp) =NaN;
    D_psd(Ns, nexp)=NaN;
    sigma_D_psd(Ns, nexp)=NaN;
end
    

%mean square displacement
maxlag_msd=300;
try
    [k_msd(Ns, nexp),sigma_k_msd(Ns, nexp),~, ~, D_msd(Ns, nexp), sigma_D_msd(Ns, nexp), ~, ~, ~, ~, gamma_msd(Ns, nexp), sigma2_gamma_msd(Ns, nexp)]=msd_nfilt(xn,T,dt*subs,maxlag_msd);
catch

    k_msd(Ns, nexp)=NaN;
    sigma_k_msd(Ns, nexp) =NaN;
    D_msd(Ns, nexp)=NaN;
    sigma_D_msd(Ns, nexp)=NaN;
end
%%%%%%%%%%%%%%%%
%especial data format for FORMA and bayesian  analysis

xr=reshape(xn, [size(xn,1)*size(xn,2),1 ]);


%FORMA

[fc_forma,D_forv,Efc_forma,sigma_D_forma(Ns, nexp)] = forma1d(xr, dt*subs);
D_forma(Ns, nexp)=mean(D_forv);
gamma_forma=kB*T./D_forma(Ns);
kfor=gamma_forma.*mean(fc_forma);

k_forma(Ns, nexp)=kfor;
sigma_k_forma(Ns, nexp)=gamma_forma.*Efc_forma;
%BAYESIAN
[k_bay(Ns, nexp), sigma_k_bay(Ns, nexp), gamma_bay(Ns, nexp), sigma_gamma_bay(Ns, nexp), D_bay(Ns, nexp), sigma_D_bay(Ns, nexp)]= bayesian(xr, dt*subs,T, a);


end

end
%%
%Average over experiments for each number of points in the trajectory


for Ns=1:length(dnn)
    %stiffness
kk_pot(Ns)=mean(k_pot(Ns, :));
sigma_kk_pot(Ns)=std(k_pot(Ns, :));

kk_eq(Ns)=mean(k_eq(Ns, :));
sigma_kk_eq(Ns)=std(k_eq(Ns, :));


kk_acf(Ns)=mean(k_acf(Ns, :));
sigma_kk_acf(Ns)=std(k_acf(Ns, :));


kk_psd(Ns)=mean(k_psd_var(Ns, :));
sigma_kk_psd(Ns)=std(k_psd_var(Ns, :));


kk_forma(Ns)=mean(k_forma(Ns, :));
sigma_kk_forma(Ns)=std(k_forma(Ns, :));


kk_bay(Ns)=mean(k_bay(Ns, :));
sigma_kk_bay(Ns)=std(k_bay(Ns, :));


kk_msd(Ns)=mean(k_msd(Ns, :));
sigma_kk_msd(Ns)=std(k_msd(Ns, :));
%%Diffusion coeficcient

DD_msd(Ns)=mean(D_msd(Ns, :));
sigma_DD_msd(Ns)=std(D_msd(Ns, :));



DD_acf(Ns)=mean(D_acf(Ns, :));
sigma_DD_acf(Ns)=std(D_acf(Ns, :));


DD_psd(Ns)=mean(D_psd(Ns, :));
sigma_DD_psd(Ns)=std(D_psd(Ns, :));


DD_forma(Ns)=mean(D_forma(Ns, :));
sigma_DD_forma(Ns)=std(D_forma(Ns, :));


DD_bay(Ns)=mean(D_bay(Ns, :));
sigma_DD_bay(Ns)=std(D_bay(Ns, :));

DD_msd(Ns)=mean(D_msd(Ns, :));
sigma_DD_msd(Ns)=std(D_msd(Ns, :));


end
save(['results_averaged_' , filename ,datestr(now, 'yyyymmddTHHMMSS'), '.mat'],'subs', 'dnn','nbins_pot' ,...
    'nw_psd','NFf', 'maxlag_msd' ,'kk_pot', 'sigma_kk_pot',...
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