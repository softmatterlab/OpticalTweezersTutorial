function [k_pot, Ek_pot,xbins, mrho, Erho, mU, EU, rho0, x_eq]=pot_nlfit(x,T,varargin)
    %POT_NLFIT   1D implementation of the POTENTIAL METHOD using non-linear fitting.
    %   POT_NLFIT(X,T, nb) generates a estimator of the stiffness k_pot
    %   for the potential method using non-linear fitting. 
    %
    %   Inputs
    %   x: time series of the position of the particle
    %   T: temperature
    %   nb: (optional) number of bins for the histogram, default is 50
    %
    %   Outputs
    %   k_pot: estimated stiffness using  non-linear fitting for the potential
    %   analysis
    %   Ek_pot: standard deviation of the stiffness
    %   xbins: x -edges for the histogram
    %   mrho: mean of the probablility distribution, average over experiments
    %   Erho: standard deviation of the probablility distribution
    %   mU: mean of the potential energy, average over experiments
    %   EU: standard deviation of the potential energy
    %   rho0: normalization factor
    %   x_eq: equilibrium position



kb=1.38064852e-23;

x = x - repmat(mean(x),size(x,1),1);

[~,Nexp]=size(x);
nb=50;

if nargin>2
    nb=varargin{1};
end

edges=linspace(min(x(:)),max(x(:)),nb);


for j=1:Nexp
    
    xx=x(:,j);

    [nc(:,j)]=histcounts(xx,edges);

end

dx=edges(2)-edges(1);

nc=nc./(sum(nc,1)*dx); % distribution

mrho=mean(nc,2);

Erho=std(nc,[],2);

Erho(Erho==0)=1;



xbins=(edges(2:end)+edges(1:end-1))/2;

% Fitting

% weights
w=1./Erho.^2;
w(isinf(w))=1;

maxbin=xbins(end); % to normalize the data to avoid "Equation is badly conditioned"

% Guess for the initial conditions for the non-linear fitting

a0=(max(mrho)*sqrt(pi))^2;

b0=0;

%Using non-linear ftting with weights

ft=fittype('sqrt(a/pi)*exp(-a*(x-b)^2)');

c=fit(xbins'/maxbin,mrho*maxbin,ft,'weights',w/maxbin^2,'StartPoint',[a0*maxbin^2,b0]);

k_pot=2*kb*T*c.a/maxbin^2; %stiffness

x_eq=c.b*maxbin^2;



cint=confint(c,0.68);  %0.68 corresponds to one standard deviation

Ek_pot=2*kb*T/maxbin^2*(cint(2,1)-cint(1,1))/2;

%estimate energy potential
logh=log(nc);

mlogh=mean(logh,2);

Elogh=std(logh,[],2);

mU=-kb*T*(mlogh-max(mlogh));

EU=-kb*T*Elogh;

rho0=sqrt(c.a/pi)*1/maxbin;  %estimation for the normalization factor
 %plots
figure(3);

clf

col1=[0,0.6,0.8];

col2=[0.9,0.4,0.2];



histogram('BinEdges',edges*1e6,'BinCounts',mrho*1e-6)

hold on

errorbar(xbins*1e6,mrho*1e-6,Erho*1e-6,'o','LineWidth',1,'Color',col2);

plot(xbins*1e6,rho0*exp(-1/2*k_pot*(xbins-x_eq).^2/(kb*T))*1e-6,'LineWidth',1,'Color',col2)

set(gca,'FontSize',16)

xlabel('$x(\mu m)$', 'Interpreter','Latex')

ylabel('$\overline{\rho^{(\textrm{ex})}}(\mu m^{-1})$','Interpreter','Latex')


title('Probability distribution, non-linear fit')
figure(4)


errorbar(xbins*1e6,mU/(kb*T),EU/(kb*T),'o','color',col1)

hold on

plot(xbins*1e6,1/2*k_pot*(xbins-x_eq).^2/(kb*T),'LineWidth',1,'color',col2)

legend({'exp','fit'},'Interpreter','latex')

set(gca,'FontSize',16)

xlabel('$x(\mu m)$','Interpreter','Latex')

ylabel('$U(k_BT)$','Interpreter','Latex')
title('Potential, non-linear fit')
%
disp('...')

disp('Potential analysis using non linear fitting')

disp(['k_pot: ' num2str(k_pot) '+-' num2str(Ek_pot)]);%agregar unidades en pN/um

