function [ k_pot, Ek_pot, xbins, mU, EU,mrho, Erho, rho0, x_eq]=pot_lfit_v1(x,T,nb)
%function [xbins, mU, EU, k_pot, Ek_pot, mhist, Ehist, h0, x0]=pot_lfit(Vx,T,nb)
% POT_LFIT   1D implementation of the POTENTIAL METHOD
%

kb=1.38064852e-23;

[~,Nexp]=size(x);

edges=linspace(min(x(:)),max(x(:)),nb);


for j=1:Nexp
    
    xx=x(:,j);

    [nc(:,j)]=histcounts(xx,edges);
end
xbins=(edges(2:end)+edges(1:end-1))/2;


dx=edges(2)-edges(1);

nc=nc./(sum(nc,1)*dx); %probability density

ind=find(nc==0);

nc(ind)=[];

xbins(ind)=[];

mrho=mean(nc,2);

Erho=std(nc,[],2);

logh=-log(nc);

mlogh=mean(logh,2);

Elogh=std(logh,[],2);

mU=kb*T*(mlogh-min(mlogh));

EU=kb*T*Elogh;

w=1./Elogh.^2;

w(isinf(w))=1;

%

max_xc=xbins(end);

xbins=xbins/max_xc; % normalization to avoid "Equation is badly conditioned"


% Approximation to a quadratic function, Using linear ftting with weights

c=fit(xbins',mlogh,'poly2','weights',w);

xbins=xbins*max_xc;

k_pot=2*kb*T*c.p1/max_xc^2;

x_eq=-c.p2/(2*c.p1)*max_xc;

rho0=exp(c.p1/max_xc^2*x_eq^2-c.p3);

cint=confint(c,0.95);

Ek_pot=2*kb*T/max_xc^2*(cint(2,1)-cint(1,1))/2;

% plots

figure(1);

clf

col1=[0,0.6,0.8];

col2=[0.9,0.4,0.2];

figure(1)
histogram('BinEdges',[xbins-dx/2,xbins(end)+dx/2]*1e6,'BinCounts',round(mrho))

hold on

errorbar(xbins*1e6,mrho,Erho,'o','LineWidth',1,'Color',col1);

plot(xbins*1e6,rho0*exp(-1/2*k_pot*(xbins-x_eq).^2/(kb*T)),'LineWidth',1,'Color',col2)

set(gca,'FontSize',16)

xlabel('$x(\mu m)$','Interpreter','Latex')

ylabel('$\langle \textrm{Counts}\rangle_{n_e}$','Interpreter','Latex')
title('Probability distribution, linear fit')
%
figure(2)

errorbar(xbins*1e6,mU/(kb*T),EU/(kb*T),'o','color',col1)

hold on

plot(xbins*1e6,1/2*k_pot*(xbins-x_eq).^2/(kb*T),'LineWidth',1,'color',col2)

legend({'exp','fit'},'Interpreter','latex')

set(gca,'FontSize',16)

xlabel('$x(\mu m)$','Interpreter','Latex')

ylabel('$U(k_BT)$','Interpreter','Latex')
title('Potential, linear fit')
%
disp('...')

disp('Potential analysis using linear fitting')

disp(['k_pot: ' num2str(k_pot) '+-' num2str(Ek_pot)]);

