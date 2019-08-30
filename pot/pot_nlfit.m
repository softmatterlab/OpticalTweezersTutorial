function [xbins, mU, EU, k_pot, Ek_pot, mhist, Ehist, h0, x0]=pot_nlfit(Vx,T,nb)
%function [xbins, mU, EU, k_pot, Ek_pot, mrho, Erho, rho0, x0]=pot_nlfit(Vx,T,nb)
% POT_NLFIT   1D implementation of the POTENTIAL METHOD using non-linear
% fitting
%

kb=1.38064852e-23;

Vx = Vx - repmat(mean(Vx),size(Vx,1),1);

[~,Nexp]=size(Vx);

edges=linspace(min(Vx(:)),max(Vx(:)),nb);


for j=1:Nexp
    
    x=Vx(:,j);

    [nc(:,j)]=histcounts(x,edges);

end

dx=edges(2)-edges(1);

nc=nc./(sum(nc,1)*dx); % distribution

mhist=mean(nc,2);

Ehist=std(nc,[],2);

Ehist(Ehist==0)=1;

logh=log(nc);

mlogh=mean(logh,2);

Elogh=std(logh,[],2);

mU=-kb*T*(mlogh-max(mlogh));

EU=-kb*T*Elogh;

xbins=(edges(2:end)+edges(1:end-1))/2;

% Fitting

% weights
w=1./Ehist.^2;

maxbin=max(abs(xbins)); % to normalize the data to avoid "Equation is badly conditioned"

% Initial conditions for the non-linear fitting

a0=(max(mhist)*sqrt(pi))^2;

% find b0

% dh=mhist-a0*sqrt(pi)*exp(-1);
% 
% ind=find(dh(1:end-1).*dh(2:end)<0);
% 
% b0=1/xbins(ind(1))^2;

b0=0;

%Using non-linear ftting with weights

ft=fittype('sqrt(a/pi)*exp(-a*(x-b)^2)');

c=fit(xbins'/maxbin,mhist*maxbin,ft,'weights',w/maxbin^2,'StartPoint',[a0*maxbin^2,b0]);

k_pot=2*kb*T*c.a/maxbin^2;

x0=c.b*maxbin^2;

h0=sqrt(c.a/pi)*1/maxbin;

cint=confint(c,0.68);

Ek_pot=2*kb*T/maxbin^2*(cint(2,1)-cint(1,1))/2;

% plots

figure(1);

clf

col1=[0,0.6,0.8];

col2=[0.9,0.4,0.2];



histogram('BinEdges',edges*1e6,'BinCounts',mhist*1e-6)

hold on

errorbar(xbins*1e6,mhist*1e-6,Ehist*1e-6,'o','LineWidth',1,'Color',col2);

plot(xbins*1e6,h0*exp(-1/2*k_pot*(xbins-x0).^2/(kb*T))*1e-6,'LineWidth',1,'Color',col2)

set(gca,'FontSize',16)

xlabel('$x(\mu m)$', 'Interpreter','Latex')

ylabel('$\overline{\rho^{(\textrm{ex})}}(\mu m^{-1})$','Interpreter','Latex')


title('Probability distribution')
figure(2)


errorbar(xbins*1e6,mU/(kb*T),EU/(kb*T),'o','color',col1)

hold on

plot(xbins*1e6,1/2*k_pot*(xbins-x0).^2/(kb*T),'LineWidth',1,'color',col2)

legend({'exp','fit'},'Interpreter','latex')

set(gca,'FontSize',16)

xlabel('$x(\mu m)$','Interpreter','Latex')

ylabel('$U(k_BT)$','Interpreter','Latex')
title('Potential')
%
disp('...')

disp('Potential analysis using non linear fitting')

disp(['k_pot: ' num2str(k_pot) '+-' num2str(Ek_pot)]);

