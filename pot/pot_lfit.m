function [xbins, mU, EU, k_pot, Ek_pot, mhist, Ehist, h0, x0]=pot_lfit(Vx,T,nb)
%function [xbins, mU, EU, k_pot, Ek_pot, mhist, Ehist, h0, x0]=pot_lfit(Vx,T,nb)
% POT_LFIT   1D implementation of the POTENTIAL METHOD
%

kb=1.38064852e-23;

[~,Nexp]=size(Vx);

edges=linspace(min(Vx(:)),max(Vx(:)),nb);


for j=1:Nexp
    
    x=Vx(:,j);

    [nc(:,j)]=histcounts(x,edges);
end

dx=edges(2)-edges(1);

%rho=nc./(sum(nc,1)*dx); %probability density

mhist=mean(nc,2);

Ehist=std(nc,[],2);

logh=-log(nc);

mlogh=mean(logh,2);

Elogh=std(logh,[],2);

ind=find(isinf(mlogh));

mlogh(ind)=[];

Elogh(ind)=[];

mhist(ind)=[];

Ehist(ind)=[];

mU=kb*T*(mlogh-min(mlogh));

EU=kb*T*Elogh;

w=1./Elogh.^2;

%
xbins=(edges(2:end)+edges(1:end-1))/2;

max_xc=max(abs(xbins));

xbins=xbins/max_xc; % normalization to avoid "Equation is badly conditioned"

xbins(ind)=[];

% Approximation to a quadratic function, Using linear ftting with weights

c=fit(xbins',mlogh,'poly2','weights',w);

xbins=xbins*max_xc;

k_pot=2*kb*T*c.p1/max_xc^2;

x0=-c.p2/(2*c.p1)*max_xc;

h0=exp(c.p1/max_xc^2*x0^2-c.p3);

cint=confint(c,0.95);

Ek_pot=2*kb*T/max_xc^2*(cint(2,1)-cint(1,1))/2;

% plots

figure(1);

clf

col1=[0,0.6,0.8];

col2=[0.9,0.4,0.2];

set(gcf,'Position',[150 300 1600 600])

axes('OuterPosition',[0 0 0.5 1])

histogram('BinEdges',[xbins-dx/2,xbins(end)+dx/2]*1e6,'BinCounts',round(mhist))

hold on

errorbar(xbins*1e6,mhist,Ehist,'o','LineWidth',1,'Color',col1);

plot(xbins*1e6,h0*exp(-1/2*k_pot*(xbins-x0).^2/(kb*T)),'LineWidth',1,'Color',col2)

set(gca,'FontSize',16)

xlabel('$x(\mu m)$','Interpreter','Latex')

ylabel('$\langle \textrm{Counts}\rangle_{n_e}$','Interpreter','Latex')

%
axes('OuterPosition',[0.5 0 0.5 1])

errorbar(xbins*1e6,mU/(kb*T),EU/(kb*T),'o','color',col1)

hold on

plot(xbins*1e6,1/2*k_pot*(xbins-x0).^2/(kb*T),'LineWidth',1,'color',col2)

legend({'exp','fit'},'Interpreter','latex')

set(gca,'FontSize',16)

xlabel('$x(\mu m)$','Interpreter','Latex')

ylabel('$U(k_BT)$','Interpreter','Latex')

%
disp('...')

disp('Potential analysis using linear fitting')

disp(['k_pot: ' num2str(k_pot) '+-' num2str(Ek_pot)]);

