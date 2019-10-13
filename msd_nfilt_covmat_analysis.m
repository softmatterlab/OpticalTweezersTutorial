% function [k_msd, Ek_msd, D_msd, ED_msd, tau, mmsd, Emsd, indc]=msd_nlfit_v1(x,T,dt,maxlag)
% %function [k_msd, Ek_msd, D_msd, ED_msd, tau, mmsd, Emsd, indc]=msd_nlfit(Vx,T,dt,maxlag)
% % MSD_NLFIT   1D implementation of the MEAN SQUARE DISPLACEMENT ANALYSIS METHOD
% % USING NON LINEAR FITTING
close all

clear;

%close all;

addpath msd



subs=1; %use a subsampled data set


    
load(['Data_positions_Fig9_1P4_S.mat']);
maxlag=1000;
dt=dt*1;
x = x - repmat(mean(x),size(x,1),1);

kb=1.38064852e-23;

[N,Nexp]=size(x);

tau=(0:maxlag-1)*dt;

for j=1:Nexp
    
    xx=x(:,j);
    
    for k=1:maxlag
    msd1(k)=mean((xx(k:end)-xx(1:end-k+1)).^2);    
    end
    
    msd(:,j)=msd1';
end

mmsd=mean(msd,2);
Y = zeros(maxlag,Nexp);
for m = 1:Nexp
    Y(:,m) = msd1(:,m) - mmsd';
end
C = Y * Y' / (maxlag - 1.0);


%%



%Emsd=std(msd,[],2)+(1e-9)^2;
Emsd=std(msd,[],2);        
%Starting points for the fitting

a0=mmsd(end); %amplitude

disp('a0')

disp(a0)

disp(Emsd(end))
%find the characteristic time

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

tau_cut=tau(1:indc);

mmsd_cut=mmsd(1:indc);

Emsd_cut=Emsd(1:indc);

max_tau=tau_cut(end);

max_mc=max(mmsd_cut);
w=inv(C);
w(isinf(w))=1;
%w=1./Emsd_cut.^2;
%ind=find(w==Inf);

%w(ind)=[];
%tau_cut(ind)=[];
%mmsd_cut(ind)=[];
% using non-linear fitting

%max_mc=1e-12;

%max_tau=1e-6;

%esta es la que estamos 



indtau=ind0;

%c=fit(tau_cut(indtau:end)'/max_tau,mmsd_cut(indtau:end)/max_mc,fittype(modelfun),'Weights',w(indtau:end)*max_mc^2,'StartPoint',[a0/max_mc,tau0/max_tau],'Display','iter','Lower',[0,0],'Upper',[inf,Inf]);%,'Robust','LAR');

%modelfun = @(b,x)(b(1)+b(2)*exp(-b(3)*x));



%opts = statset('nlinfit');
%opts.RobustWgtFun = 'bisquare';
guess = [a0/max_mc,tau0/max_tau];
%[beta,R,J,CovB,MSE] = nlinfit(tau_cut(indtau:end)'/max_tau,mmsd_cut(indtau:end)/max_mc,modelfun,beta0);%, 'Weights',w(indtau:end)*max_mc^2);%,opts);

%%
chi2 = @(p, t, y,C ) (y - f(t, p))' * inv(C') * (y - f(t, p));

% Do the parameter fitting
[params, chi2_min] = fminunc(@(par)(chi2(par,tau_cut(indtau:length(tau_cut))'/max_tau , mmsd_cut(indtau:length(mmsd_cut))/max_mc,   C)), guess);

%cint = nlparci(beta,R,'covar',w)    
%xrange = min(x):1e-9:max(x);
%[ypred,delta, kk] = nlpredci(modelfun,xrange,beta,R,'Covar',CovB,...
                       %  'MSE',MSE,'SimOpt','on');
%%
tau0=beta(2)*max_tau;

a0=beta(1)*max_mc;

k_msd=2*kb*T/a0;


gamma_msd=k_msd*tau0;

D_msd=kb*T/gamma_msd;

%cint_uu=confint(beta,0.95);

%Ek_msd=kb*T/a0^2*(cint(2,1)-cint(1,1))/2*max_mc;
Ek_msd=abs(kb*T/(a0)^2*(cint(1,2)-cint(1,1)))*max_mc;
%Ek_msd_uu=kb*T/(a0)^2*(cint_uu(2,1)-cint_uu(1,1))*max_mc


Etau=(cint(2,2)-cint(1,2))/2*max_tau;

ED_msd=kb*T/(k_msd^2*tau0)*Ek_msd+kb*T/(k_msd*tau0^2)*Etau;

Egamma_msd=Ek_msd*tau0+k_msd*Etau;
% plot
figure(1)

clf

set(gcf,'Position',[150 300 800 600])

axes('OuterPosition',[0 0 1 1])

errorbar(tau(1:1:end),mmsd(1:1:end)*1e12,Emsd(1:1:end)*1e12,'ob','LineWidth',1)

hold on

plot(tau_cut(1:1:end),a0*(1-exp(-tau_cut(1:1:end)/tau0))*1e12,'b')

plot([ntaus*tau0,ntaus*tau0],[0.8*min(mmsd),1.3*max(mmsd)]*1e12,'--k')

axis([tau(1),tau(end),0.8*min(mmsd)*1e12,1.2*max(mmsd)*1e12])

text(ntaus*tau0,0.05*max(mmsd)*1e12,[num2str(ntaus),'$\tau_0$'],'Interpreter','latex','FontSize',16)

set(gca,'FontSize',16)

xlabel('$\tau$','Interpreter','latex')

ylabel('$MSD_x(\mu \textrm{m}^2)$','Interpreter','latex')

%   
disp('...')

disp('MSD by non-linear fitting')

disp(['tau0_msd: ' num2str(tau0) '+-' num2str(Etau)])

disp(['k_msd: ' num2str(k_msd*1e6) '+-' num2str(Ek_msd*1e6)])

disp(['D_msd: ' num2str(D_msd) '+-' num2str(ED_msd)])

disp(['gamma_msd: ' num2str(gamma_msd) '+-' num2str(Egamma_msd)])



