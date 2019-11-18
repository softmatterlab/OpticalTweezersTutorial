function [k_acf, sigma_k_acf, D_acf, sigma_D_acf,gamma_acf, sigma_gamma_acf,tau, mc, Ec,indc, tau0_exp,sigma_tau0_exp, c0_exp]=acf_lfit(Vx,T,dt)
%function [k_acf, Ek_acf, D_acf, ED_acf, tau, mc, Ec]=acf_lfit(Vx,T,dt)
% ACF_LFIT   1D implementation of the AUTOCORRELATION ANALYSIS METHOD
% USING LINEAR FITTING
addpath wlsice
 
Vx = Vx - repmat(mean(Vx),size(Vx,1),1);
 
kb=1.38064852e-23;
 
[N,Nexp]=size(Vx);
 
tau=(0:N-1)*dt;
acf=zeros(Nexp, N);
for j=1:Nexp 
    x=Vx(:,j);
    c=xcorr(x,'Unbiased');
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
 
ntaus=2.5;
 
indc=round(ntaus*ind); % consider only ntaus times the characteristic time in the fitting
acf_cut=acf(:, 1:3:indc);
tau_cut=tau(1:3:indc);
max_tau=max(tau_cut);
 
% using non-linear fitting
max_mc=max(max(acf_cut));
 
guess=[1,1];
[params, sigma, chi2_min, C] = wlsice(tau_cut/tau0, log(abs(acf_cut/max_mc)), guess, 1);

tau0_exp=tau0*params(1);
sigma_tau0_exp=tau0*sigma(1);

c0_exp=exp(params(2))*abs(max_mc);
sigma_c0_exp=exp(params(2))*abs(max_mc)*sigma(2);

k_acf=kb*T/c0_exp;
sigma_k_acf=(kb*T/c0_exp^2)*sigma_c0_exp;

gamma_acf=k_acf*tau0_exp;
sigma_gamma_acf=k_acf*sigma_tau0_exp+tau0_exp*sigma_k_acf;

D_acf=kb*T/gamma_acf;

sigma_D_acf=(kb*T/(gamma_acf^2))*sigma_gamma_acf;


end

