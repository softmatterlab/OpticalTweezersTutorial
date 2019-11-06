% USING LINEAR FITTING

close all;

clear all;

%close all;

addpath acf
addpath wlsice

    
load('Data_positions_Fig9_1P6_S.mat');



kB=1.38064852e-23;
addpath wlsice
 
x = x - repmat(mean(x),size(x,1),1);
 
kb=1.38064852e-23;
 
[N,Nexp]=size(x);
 
tau=(0:N-1)*dt;
acf=zeros(Nexp, N);
for j=1:Nexp 
    xpos=x(:,j);
    c=xcorr(xpos,'Unbiased');
    c=c(N:end);
    acf(j,:)=c;
end
 
mc=mean(acf,1);
 
Ec=std(acf,[],1);
 
% first approximation to define the starting points and the significative
% points in the fitting
 
c0=mc(1); %amplitude
 
ctau=c0*exp(-1); 
 
dc=mc-ctau;
 
%find the characteristic time
 
ind=find(dc(1:end-1).*dc(2:end)<0);
 
tau0=tau(ind(1));
 
ntaus=4;
 
indc=round(ntaus*ind); % consider only ntaus times the characteristic time in the fitting
acf_cut=acf(:, 1:indc);
tau_cut=tau(1:indc);
max_tau=max(tau_cut);
 
% using non-linear fitting
max_mc=max(max(acf_cut));
 
guess=[1,1];
[params, sigma, chi2_min, C] = wlsice(tau_cut/tau0, log(abs(acf_cut/max_mc)), guess, 1);
sigma 
params
tau0_exp=tau0*params(1)
sigma_tau0_exp=tau0*sigma(1)

c0_exp=exp(params(2))*abs(max_mc)
sigma_c0_exp=exp(params(2))*abs(max_mc)*sigma(2)
k_acf=kb*T/c0_exp

sigma_k_acf=(kb*T/c0_exp^2)*sigma_c0_exp
%1/exp(params(2))*sigma(2)/max_mc
%Ek_acf=kb*T/exp(params(2))*sigma(2)/max_mc;
gamma_acf=k_acf*tau0_exp;
sigma_gamma_acf=k_acf*sigma_tau0_exp+tau0_exp*sigma_k_acf;

D_acf=kb*T/gamma_acf;

sigma_D_acf=(kb*T/(gamma_acf^2))*sigma_gamma_acf;
%Ek_acf=kb*T*exp(-params(2)*log(c0))*log(c0)*(sigma(2));Isaac


%ED_acf=kb*T/(k_acf^2*tau0)*Ek_acf+kb*T/(k_acf*tau0^2)*abs(1/cint(2,1)-1/cint(1,1))/2;
ED_acf=log(abs(max_mc))*max_mc/tau0*(exp(params(2))*sigma(2)*params(1)+sigma(1)*exp(params(2)));
gamma_exp=kb*T/D_acf;
 
sigma2_gamma_exp=kb*T/(D_acf^2)*ED_acf;
% % plot
% figure(1)
% 
% 
% clf
% 
% set(gcf,'Position',[150 300 1600 600])
% 
% axes('OuterPosition',[0 0 1 1])
% 
%
%errorbar(tau(1:20:6*indc),mc(1:20:6*indc)*1e12,Ec(1:20:6*indc)*1e12,'ob','LineWidth',1)
% 
 %hold on
% 
 %plot(tau(1:20:6*indc),c0_exp*exp(tau(1:20:6*indc)/tau0_exp)*1e12,'b')
% 
% xlabel('$\tau$','Interpreter','latex')
% 
% ylabel('$C_x(\mu \textrm{m}^2)$','Interpreter','latex')
% 
% %
% disp('...')
% % 
%  disp('Autocorrelation function analysis by linear fitting')
% % 
  disp(['k_acf: ' num2str(k_acf*1e6) '+-' num2str(sigma_k_acf*1e6)])
% % 
disp(['D_acf: ' num2str(D_acf*1e12) '+-' num2str(sigma_D_acf*1e12)])
% 
 disp(['gamma_acf:' num2str(gamma_acf*1e9) '+-' num2str(sigma_gamma_acf*1e9)])

