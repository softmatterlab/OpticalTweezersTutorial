function [k_acf, Ek_acf, D_acf, ED_acf, tau, mc, Ec, indc]=acf_nlfit(Vx,T,dt)
%function [k_acf, Ek_acf, D_acf, ED_acf, tau, mc, Ec]=acf_nlfit(Vx,T,dt)
% ACF_NLFIT   1D implementation of the AUTOCORRELATION ANALYSIS METHOD
% USING NON LINEAR FITTING

Vx = Vx - repmat(mean(Vx),size(Vx,1),1);

kb=1.38064852e-23;

[N,Nexp]=size(Vx);

tau=(0:N-1)*dt;

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

ntaus=2;

indc=ntaus*ind; % consider only ntaus times the characteristic time in the fitting

tau_cut=tau(1:indc);

mc_cut=mc(1:indc);

Ec_cut=Ec(1:indc);
acf_cut=acf(:, 1:indc);

max_tau=max(tau_cut);

max_mc=max(mc_cut);

w=1./Ec_cut.^2;

% using non-linear fitting

%ft=fittype('a*exp(-x/b)');
guess=[c0/max_mc,tau0/max_tau];
[params, sigma, chi2_min, C] = wlsice(tau_cut/max_tau, acf_cut/max_mc, guess);


%c=fit(tau_cut'/max_tau,mc_cut/max_mc,ft,'Weights',w*max_mc,'StartPoint',[c0/max_mc,tau0/max_tau]);

tau0=params(2)*max_tau

c0=params(1)*max_mc;

%cint=[params(1)-sigma(1)/2 params(2)-sigma(2)/2; params(1)+sigma(1)/2 params(2)+sigma(2)/2];


k_acf=kb*T/c0

D_acf=kb*T/(k_acf*tau0);

%cint=confint(c,0.95);

Ek_acf=kb*T/c0^2*(sigma(2))/2*max_mc

ED_acf=kb*T/(k_acf^2*tau0)*Ek_acf+kb*T/(k_acf*tau0^2)*(sigma(2))/2*max_tau;

% plot
figure(1)
clf

errorbar(tau(1:20:6*indc),mc(1:20:6*indc)*1e12,Ec(1:20:6*indc)*1e12,'ob','LineWidth',1)

hold on

plot(tau(1:20:6*indc),c0*exp(-tau(1:20:6*indc)/tau0)*1e12,'b')

xlabel('$\tau$','Interpreter','latex')

ylabel('$C_x(\mu \textrm{m}^2)$','Interpreter','latex')

%   
disp('...')

disp('Autocorrelation function analysis by non-linear fitting')

disp(['k_acf: ' num2str(k_acf*1e6) '+-' num2str(Ek_acf*1e6)])

disp(['D_acf: ' num2str(D_acf) '+-' num2str(ED_acf)])

disp(['gamma_acf:' num2str(kb*T/D_acf) '+-' num2str(kb*T/D_acf^2*ED_acf)])
disp(['tau_0:' num2str(tau0)]);

