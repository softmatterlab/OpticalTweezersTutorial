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
Nexp=5;
subs=1;
dnn=floor((1e2:5*5e2:1e6)/Nexp);
disp(size(dnn));
%dnn=1e6/Nexp;
disp(length(dnn));
disp(dnn(end))
%%
for nne=0:Nexp-1
disp(nne);

    %first we create the Nexp different experiments
xs=x(1+nne*maxn:1+(nne+1)*maxn, :);

disp(size(xs));


for Ns=1:length(dnn)
%then we start to change the number of points in the trajectory
disp(Ns);


xn =xs(1:dnn(Ns), :);
N=size(xn);
disp(size(xn));  
%Potential
nbins_pot=50;
try

    [k_pot(Ns, nne+1), sigma_k_pot(Ns, nne+1), ~, ~, ~, ~, ~, ~, ~,~]=pot_lfit(xn,T,nbins_pot);
catch
    k_pot(Ns, nne+1)=NaN;
    sigma_k_pot(Ns, nne+1) =NaN;
end

%Equipartition
try
    [k_eq(Ns, nne+1), sigma_k_eq(Ns, nne+1)]=eq1d(xn,T);
catch
    k_eq(Ns, nne+1)=NaN;
    sigma_k_eq(Ns, nne+1) =NaN;
end
%autocorrelation function
try
    
    [k_acf(Ns, nne+1), sigma_k_acf(Ns, nne+1), D_acf(Ns, nne+1), sigma_D_acf(Ns, nne+1),...
    gamma_acf(Ns, nne+1), sigma_gamma_acf(Ns, nne+1),~, ~, ~,~, ~, ~]=acf_nlfit(xn,T,dt*subs);
catch
    k_acf(Ns, nne+1)=NaN;
    sigma_k_acf(Ns, nne+1) =NaN;
    D_acf(Ns, nne+1)=NaN;
    sigma_D_acf(Ns, nne+1)=NaN;
    gamma_acf(Ns, nne+1)=NaN;
    sigma_gamma_acf(Ns, nne+1)=NaN;
end
%power spectral density
% number of windows
nw_psd=500;
nw=round(1*size(xn,1)/nw_psd);
NFf=1/4;
try
    

[fc_psd,D_psd(Ns, nne+1),sigma_fc_psd,sigma_D_psd(Ns, nne+1),f,XX,fw_mean,Pk,sigma_Pk,fcut]=psd_lfit(xn,dt*subs,nw,NFf);

gamma_psd(Ns, nne+1)=kB*T/D_psd(Ns, nne+1);

%sigma_gamma_psd(Ns, nne+1)=kb*T./D_psd^2sigma_D_psd(Ns, nne+1);



% estimation of k using the estimated gamma
k_psd(Ns, nne+1)=2*pi*gamma_psd(Ns, nne+1)*fc_psd;

%sigma_k_psd=2*pi*gamma_psd(Ns, nne+1).*sigma_fc_psd+2*pi*fc_psd*sigma_gamma_psd;
%     [mfc_psd(Ns, nne+1),D_psd(Ns, nne+1),Efc_psd,sigma_D_psd(Ns, nne+1),f,XX,fw_mean,Pk,EPk,fcut]=psd_lfit(xn,dt*subs,nw,NFf);
%     
%     gamma_psd(Ns, nne+1)=kB*T./mean(D_psd(Ns, nne+1));
%     
%     sigma_gamma_psd(Ns, nne+1)=kB*T./mean(D_psd(Ns, nne+1))^2*mean(sigma_D_psd(Ns, nne+1));
%     
%     k_psd(Ns, nne+1)=2*pi*gamma*mean(mfc_psd(Ns, nne+1));
%     
%     sigma_k_psd=2*pi*gamma*Efc_psd;
%     
%     % estimation of k using the estimated gamma
%     k_psd_var(Ns, nne+1)=2*pi*gamma_psd(Ns, nne+1).*mean(mfc_psd(Ns, nne+1))
%     
%     sigma_k_psd_var(Ns, nne+1)=2*pi*gamma_psd.*Efc_psd+2*pi*mean(mfc_psd(Ns, nne+1))*sigma_gamma_psd(Ns, nne+1);
catch ME
    %rethrow(ME)
    k_psd(Ns, nne+1)=NaN;
    sigma_k_psd(Ns, nne+1) =NaN;
    D_psd(Ns, nne+1)=NaN;
    sigma_D_psd(Ns, nne+1)=NaN;
    gamma_psd(Ns, nne+1)=NaN;
    sigma_gamma_psd(Ns, nne+1)=NaN;
end
    

%mean square displacement
maxlag_msd=1000;
try
    [k_msd(Ns, nne+1),sigma_k_msd(Ns, nne+1),~, ~, D_msd(Ns, nne+1), sigma_D_msd(Ns, nne+1), ~, ~, ~, ~, gamma_msd(Ns, nne+1), sigma2_gamma_msd(Ns, nne+1)]=msd_nfilt(xn,T,dt*subs,maxlag_msd);
catch

    k_msd(Ns, nne+1)=NaN;
    sigma_k_msd(Ns, nne+1) =NaN;
    D_msd(Ns, nne+1)=NaN;
    sigma_D_msd(Ns, nne+1)=NaN;
    gamma_msd(Ns, nne+1)=NaN;
    sigma2_gamma_msd(Ns, nne+1)=NaN;
end
%%%%%%%%%%%%%%%%
%especial data format for FORMA and bayesian  analysis

xr=reshape(xn, [size(xn,1)*size(xn,2),1 ]);
nsubs=3;

%FORMA

[fc_forma,D_forv,Efc_forma(Ns, nne+1),sigma_D_forma(Ns, nne+1)] = forma1d(xr, dt, nsubs);
D_forma(Ns, nne+1)=mean(D_forv);
gamma_forma(Ns, nne+1)=kB*T./D_forma(Ns);
kfor=gamma_forma(Ns, nne+1).*mean(fc_forma);

k_forma(Ns, nne+1)=kfor;
sigma_k_forma(Ns, nne+1)=gamma_forma(Ns, nne+1).*Efc_forma(Ns, nne+1);
%BAYESIAN
[k_bay(Ns, nne+1), sigma_k_bay(Ns, nne+1), gamma_bay(Ns, nne+1), sigma_gamma_bay(Ns, nne+1), D_bay(Ns, nne+1), sigma_D_bay(Ns, nne+1)]= bayesian(xr, dt,T, a, nsubs);


end

end
%%
%Average over experiments for each number of points in the trajectory
kk_pot=mean(k_pot, 2);
sigma_kk_pot=std(k_pot,[], 2);

kk_eq=mean(k_eq, 2);
sigma_kk_eq=std(k_eq,[], 2);

kk_acf=mean(k_acf, 2);
 sigma_kk_acf=std(k_acf,[], 2);

 
kk_psd=mean(k_psd, 2);
sigma_kk_psd=std(k_psd,[], 2);


kk_forma=mean(k_forma, 2);
sigma_kk_forma=std(k_forma,[], 2);


kk_msd=mean(k_msd, 2);
sigma_kk_msd=std(k_msd,[], 2);


% kk_bay(Ns)=mean(k_bay(Ns, :));
% sigma_kk_bay(Ns)=std(k_bay(Ns, :));

% for Ns=1:length(dnn)
%     %stiffness
% kk_pot(Ns)=mean(k_pot(Ns, :));
% sigma_kk_pot(Ns)=std(k_pot(Ns, :));
% 
% kk_eq(Ns)=mean(k_eq(Ns, :));
% sigma_kk_eq(Ns)=std(k_eq(Ns, :));
% 
% 
% kk_acf(Ns)=mean(k_acf(Ns, :));
% sigma_kk_acf(Ns)=std(k_acf(Ns, :));
% 
% 
% kk_psd(Ns)=mean(k_psd(Ns, :));
% sigma_kk_psd(Ns)=std(k_psd(Ns, :));
% 
% 
% kk_forma(Ns)=mean(k_forma(Ns, :));
% sigma_kk_forma(Ns)=std(k_forma(Ns, :));
% 
% 
% kk_bay(Ns)=mean(k_bay(Ns, :));
% sigma_kk_bay(Ns)=std(k_bay(Ns, :));
% 
% 
% kk_msd(Ns)=mean(k_msd(Ns, :));
% sigma_kk_msd(Ns)=std(k_msd(Ns, :));
% %%Diffusion coeficcient
% 
% DD_msd(Ns)=mean(D_msd(Ns, :));
% sigma_DD_msd(Ns)=std(D_msd(Ns, :));
% 
% 
% 
% DD_acf(Ns)=mean(D_acf(Ns, :));
% sigma_DD_acf(Ns)=std(D_acf(Ns, :));
% 
% 
% DD_psd(Ns)=mean(D_psd(Ns, :));
% sigma_DD_psd(Ns)=std(D_psd(Ns, :));
% 
% 
% DD_forma(Ns)=mean(D_forma(Ns, :));
% sigma_DD_forma(Ns)=std(D_forma(Ns, :));
% 
% 
% DD_bay(Ns)=mean(D_bay(Ns, :));
% sigma_DD_bay(Ns)=std(D_bay(Ns, :));
% 
% DD_msd(Ns)=mean(D_msd(Ns, :));
% sigma_DD_msd(Ns)=std(D_msd(Ns, :));
% 
% 
% 
% %Friction coefficent
% 
% %%Diffusion coeficcient
% 
% gg_msd(Ns)=mean(gamma_msd(Ns, :));
% sigma_gg_msd(Ns)=std(gamma_msd(Ns, :));
% 
% 
% 
% gg_acf(Ns)=mean(gamma_acf(Ns, :));
% sigma_gg_acf(Ns)=std(gamma_acf(Ns, :));
% 
% 
% gg_psd(Ns)=mean(gamma_psd(Ns, :));
% sigma_gg_psd(Ns)=std(gamma_psd(Ns, :));
% 
% 
% gg_forma(Ns)=mean(gamma_forma(Ns, :));
% sigma_gg_forma(Ns)=std(gamma_forma(Ns, :));
% 
% 
% gg_bay(Ns)=mean(gamma_bay(Ns, :));
% sigma_gg_bay(Ns)=std(gamma_bay(Ns, :));



save(['results_averaged_' , filename ,datestr(now, 'yyyymmddTHHMMSS'), '.mat'],'subs', 'dnn','nbins_pot' ,...
    'nw_psd','NFf', 'maxlag_msd' ,'kk_pot', 'sigma_kk_pot',...
    'kk_eq', 'sigma_kk_eq','kk_acf', 'sigma_kk_acf', 'kk_psd', 'sigma_kk_psd',...
    'kk_msd', 'sigma_kk_msd', 'kk_forma', 'sigma_kk_forma', 'kk_bay', 'sigma_kk_bay', ...
    'DD_acf', 'sigma_DD_acf', 'DD_psd', 'sigma_DD_psd', 'DD_msd', 'sigma_DD_msd','DD_forma',...
    'sigma_DD_forma', 'DD_bay', 'sigma_DD_bay','gg_msd', 'sigma_gg_msd', 'gg_acf' , 'sigma_gg_acf',....
    'gg_psd', 'sigma_gg_psd', 'gg_forma', 'sigma_gg_forma', 'gg_bay')
%%

%dnn=dnn(1:50);
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


figure(5)
%load('results_comparison.mat')

semilogx(dnn, gg_acf, 'DisplayName', 'ACF', 'Color', 'cyan', 'LineWidth', 3)
hold on
semilogx(dnn, gg_psd, 'DisplayName', 'PSD', 'Color', 'yellow', 'LineWidth', 3)
semilogx(dnn, gg_msd, 'DisplayName', 'MSD', 'Color', 'black', 'LineWidth', 3)

semilogx(dnn, gg_forma, 'DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 3)
hold on
semilogx(dnn, gg_bay, 'DisplayName', 'BAYESIAN', 'Color', [0.1 0.5 0.7], 'LineWidth', 3)
legend
xlabel('$N_s$', 'Interpreter', 'Latex', 'FontSize',30)
ylabel('$k( p \rm N  \mu \rm m^{-1})$', 'Interpreter', 'Latex', 'FontSize',30)
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
hold off