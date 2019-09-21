function [k_pot, Ek_pot, xbins, mrho, Erho, mU, EU, rho0, x_eq]=pot_nlfit(x,T,varargin)
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

%define the histogram

for j=1:Nexp
    
    xx=x(:,j);

    [nc(:,j)]=histcounts(xx,edges);

end

dx=edges(2)-edges(1);

xbins=(edges(2:end)+edges(1:end-1))/2;

nc=nc./(sum(nc,1)*dx); % probability density

mrho=mean(nc,2);

[ind]=find(mrho==0);

xbins(ind)=[];

Erho=std(nc,[],2);

Erho(Erho==0)=1;

% weights

w=1./Erho.^2;

w(isinf(w))=1;

maxbin=xbins(end); 

xbins=xbins/maxbin; % normalization to avoid "Equation is badly conditioned"

% Guess for the initial conditions for the non-linear fitting

a0=(max(mrho)*sqrt(pi))^2;

b0=0;
ind=find(mrho==Inf);

mrho(ind)=[];

xbins(ind)=[];

%Using non-linear ftting with weights

ft=fittype('sqrt(a/pi)*exp(-a*(x-b)^2)');

c=fit(xbins',mrho*maxbin,ft,'weights',w/maxbin^2,'StartPoint',[a0*maxbin^2,b0]);

%return to original variables after fit

xbins=xbins*maxbin;

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


disp('...')

disp('Potential analysis using non linear fitting')

disp(['k_pot: ' num2str(k_pot) '+-' num2str(Ek_pot)]);%agregar unidades en pN/um

