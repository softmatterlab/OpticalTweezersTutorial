function [k_msd,sigma_k_msd, tau0, sigma2_tau0, D_msd, ED_msd, tau, mmsd, Emsd, indc, gamma_msd, sigma2_gamma_msd]=msd_nfilt(x,T,dt,maxlag)
%function [k_msd, Ek_msd, D_msd, ED_msd, tau, mmsd, Emsd, indc]=msd_nlfit(Vx,T,dt,maxlag)
% MSD_NLFIT   1D implementation of the MEAN SQUARE DISPLACEMENT ANALYSIS METHOD
% USING NON LINEAR FITTING
addpath wlsice
x = x - repmat(mean(x),size(x,1),1);

kb=1.38064852e-23;

[~,Nexp]=size(x);

tau=(0:maxlag-1)*dt;

for j=1:Nexp
    
    xx=x(:,j);
    
    for k=1:maxlag
    msd1(k)=mean((xx(k:end)-xx(1:end-k+1)).^2);    
    end
    
    msd(j,:)=msd1;
end

mmsd=mean(msd,1);

%Emsd=std(msd,[],2)+(1e-9)^2;
Emsd=std(msd,[],1);        
%Starting points for the fitting

a0=mmsd(end); %amplitude




%find the characteristic time

atau=a0*(1-exp(-1)); 

dc=mmsd-atau;

ind=find(dc(1:end-1).*dc(2:end)<0);

indtau=ind(1);

tau0=tau(indtau);


ntaus=6;

indc=ntaus*ind; % consider only ntaus times the characteristic time in the fitting

if indc>maxlag
    indc=maxlag;
end

ind0=floor(ind/3);

tau_cut=tau(1:indc);

mmsd_cut=mmsd(1:indc);
msd_cut_exp=msd(:, 1:indc);

Emsd_cut=Emsd(1:indc);

max_tau=tau_cut(end);

max_mc=max(mmsd_cut);

%w=1./Emsd_cut.^2;
%ind=find(w==Inf);

%w(ind)=[];
%tau_cut(ind)=[];
%mmsd_cut(ind)=[];
% using non-linear fitting

%max_mc=1e-12;

%max_tau=1e-6;

%ftmodelfun=@(b,x)(b(1))*(1-exp(-x/b(2)));

indtau=ind0;

%c=fit(tau_cut(indtau:end)'/max_tau,mmsd_cut(indtau:end)/max_mc,fittype(ft),'Weights',w(indtau:end)*max_mc^2,'StartPoint',[a0/max_mc,tau0/max_tau],'Display','iter','Lower',[0,0],'Upper',[inf,Inf]);%,'Robust','LAR');

guess=[a0/max_mc,tau0/max_tau];
%[params, sigma, chi2_min, C] = wlsice(time, trajectories, guess, opt)
[params, sigma, chi2_min, C] = wlsice(tau_cut(indtau:end)/max_tau, msd_cut_exp(:, indtau:end)/max_mc,guess, 3);
%[params, chi2_min, C, sigma] = fit_cov(tau_cut(indtau:end)'/max_tau, mmsd_cut(indtau:end)/max_mc,guess);
%%

tau0=params(2)*max_tau;
sigma2_tau0=sigma(2)/2*max_tau;
a0=params(1)*max_mc;

k_msd=2*kb*T/a0;


gamma_msd=k_msd*tau0;

D_msd=kb*T/gamma_msd;

%cint=confint(c,0.95);

%Ek_msd=kb*T/a0^2*(cint(2,1)-cint(1,1))/2*max_mc;
sigma_k_msd=kb*T/(a0)^2*(sigma(1))*max_mc;

Etau=(sigma(2))/2*max_tau;

ED_msd=kb*T/(k_msd^2*tau0)*sigma_k_msd+kb*T/(k_msd*tau0^2)*Etau;

sigma2_gamma_msd=sigma_k_msd*tau0+k_msd*Etau;
% disp('k_msd')
% disp([num2str(k_msd*1e6) '+-'  num2str(sigma_k_msd*1e6)])
% disp('tau0')
% disp([num2str(tau0) '+-'  num2str(sigma2_tau0)])
% disp('gamma')
% disp([num2str(gamma_msd) '+-'  num2str(sigma2_gamma_msd)])

