% Licensed under the terms of the GPL v3. See LICENSE for details.
clear all 
close all
addpath msd

subs=1; %use a subsampled data set


load(['Data_positions_Fig9_1P4_S.mat']);
xx=x(1:subs:size(x,1),:);

xx = xx - repmat(mean(xx),size(xx,1),1);
maxlag=500;
dt=dt*subs;

kb=1.38064852e-23;

[N,Nexp]=size(xx);
tau=(0:maxlag-1)*dt;
max(tau)
for j=1:Nexp
    
    xxe=xx(:,j);
    
    for k=1:maxlag
    msd1(k)=mean((xxe(k:end)-xxe(1:end-k+1)).^2);    
    end
    
    msd(:,j)=msd1';
end

mmsd=mean(msd,2);

Emsd=std(msd,[],2);  

errorbar(tau,mmsd, Emsd)
hold on
%%
a0=mmsd(end); %amplitude

disp('a0')

disp(a0)

disp(Emsd(end))


atau=a0*(1-exp(-1)); 

dc=mmsd-atau;

ind=find(dc(1:end-1).*dc(2:end)<0);

indtau=ind(1);

tau0=tau(indtau);
disp('tau0')
disp(tau0)

ntaus=6;

indc=ntaus*ind; % consider only ntaus times the characteristic time in the fitting

if indc>maxlag
    indc=maxlag;
end

ind0=floor(ind/4);

tau_cut=tau(ind0:indc);

mmsd_cut=mmsd(ind0:indc);

Emsd_cut=Emsd(ind0:indc);

max_tau=tau_cut(end);

max_mc=max(mmsd_cut);


guess = [a0/max_mc,tau0/max_tau];


[params, chi2_min, C, sigma] = fit_cov(tau_cut/max_tau, mmsd_cut/max_mc,guess);
%%

%%



tau0=params(2)*max_tau;

a0=params(1)*max_mc;

k_msd=2*kb*T/a0;


gamma_msd=k_msd*tau0;

D_msd=kb*T/gamma_msd;


Ek_msd=kb*T/(a0)^2*(sigma(1))*max_mc;
Etau=(sigma(2))/2*max_tau;

ED_msd=kb*T/(k_msd^2*tau0)*Ek_msd+kb*T/(k_msd*tau0^2)*Etau;

Egamma_msd=Ek_msd*tau0+k_msd*Etau;

msd_mod=a0*(1-exp(-tau/tau0));
plot(tau, msd_mod)

disp('...')

disp('MSD by non-linear fitting')

disp(['tau0_msd: ' num2str(tau0) '+-' num2str(Etau)])

disp(['k_msd: ' num2str(k_msd*1e6) '+-' num2str(Ek_msd*1e6)])

disp(['D_msd: ' num2str(D_msd*1e12) '+-' num2str(ED_msd*1e12)])

disp(['gamma_msd: ' num2str(gamma_msd) '+-' num2str(Egamma_msd)])