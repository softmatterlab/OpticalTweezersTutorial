function [k_acf, Ek_acf, D_acf, ED_acf, tau, mc, Ec,indc, tau0_exp, c0_exp, acf, C]=acf_lfit(Vx,T,dt)
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
disp(size(acf))
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

ntaus=1.5;

indc=round(ntaus*ind); % consider only ntaus times the characteristic time in the fitting
acf_cut=acf(:, 1:3:indc);
tau_cut=tau(1:3:indc);
max_tau=max(tau_cut);

% using non-linear fitting
max_mc=max(max(acf_cut))



guess=[1,1];
[params, sigma, chi2_min, C] = wlsice(tau_cut/tau0, log(abs(acf_cut))/log(abs(max_mc)), guess, 'acf_lf');
params
sigma

tau0_exp=tau0/params(1)/log(abs(max_mc))

c0_exp=exp(params(2)*log(max_mc));

k_acf=kb*T/c0_exp;

D_acf=kb*T/(k_acf*tau0);


%Ek_acf=kb*T*exp(-params(2)*log(c0))*log(c0)*(sigma(2));
Ek_acf=kb*T/exp(params(2))*sigma(2)/max_mc;

%ED_acf=kb*T/(k_acf^2*tau0)*Ek_acf+kb*T/(k_acf*tau0^2)*abs(1/cint(2,1)-1/cint(1,1))/2;
ED_acf=0;
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
errorbar(tau(1:20:6*indc),mc(1:20:6*indc)*1e12,Ec(1:20:6*indc)*1e12,'ob','LineWidth',1)
% 
 hold on
% 
 plot(tau(1:20:6*indc),c0_exp*exp(tau(1:20:6*indc)/tau0_exp)*1e12,'b')
% 
% xlabel('$\tau$','Interpreter','latex')
% 
% ylabel('$C_x(\mu \textrm{m}^2)$','Interpreter','latex')
% 
% %
disp('...')
% 
 disp('Autocorrelation function analysis by linear fitting')
% 
 disp(['k_acf: ' num2str(k_acf*1e6) '+-' num2str(Ek_acf*1e6)])
% 
% disp(['D_acf: ' num2str(D_acf) '+-' num2str(ED_acf)])
% 
% disp(['gamma_acf:' num2str(kb*T/D_acf) '+-' num2str(kb*T/D_acf^2*ED_acf)])
end
