function [k_pot, sigma2_k_pot, x_alpha, mrho, sigma2_rho, mU, sigma2_U, rho0, x_eq]=pot_lfit(x,T,varargin)
    %POT_LFIT   1D implementation of the POTENTIAL METHOD using linear fitting.
    %   POT_LFIT(X,T, nb) generates a estimator of the stiffness k_pot
    %   for the potential method using linear fitting. 
    %
    %   Inputs
    %   x: time series of the position of the particle
    %   T: temperature
    %   P: (optional) number of bins for the histogram, default is 50
    %
    %   Outputs
    %   k_pot: estimated stiffness using  non-linear fitting for the potential
    %   analysis
    %   sigma2_k_pot: standard deviation squared of the stiffness
    %   x_alpha: x -edges for the histogram
    %   mrho: mean of the probablility distribution, average over experiments
    %   sigma2_rho: standard deviation squared of the probablility distribution
    %   mU: mean of the potential energy, average over experiments
    %   sigma2_U: standard deviation squared of the potential energy
    %   rho0: normalization factor
    %   x_eq: equilibrium position

kb=1.38064852e-23;

x = x - repmat(mean(x),size(x,1),1);

[~,Nexp]=size(x);

P=50;

if nargin>2
    P=varargin{1};
end

edges=linspace(min(x(:)),max(x(:)),P);

%define the histogram

for j=1:Nexp
    
    xx=x(:,j);

    [frequency(:,j)]=histcounts(xx,edges);
    
end

dx=edges(2)-edges(1);

x_alpha=(edges(2:end)+edges(1:end-1))/2;

frequency=frequency./(sum(frequency,1)*dx); % after normalization ot becomes the probability density

mrho=mean(frequency,2);

[ind]=find(mrho==0);

x_alpha(ind)=[];

sigma2_rho=std(frequency,[],2);

sigma2_rho(sigma2_rho==0)=1;

% weights for fitting

w=1./sigma2_rho.^2;

w(isinf(w))=1;

maxbin=x_alpha(end);

x_alpha=x_alpha/maxbin; % normalization to avoid "Equation is badly conditioned"

logh=-log(frequency);

mlogh=mean(logh,2);

[ind]=find(mlogh==Inf);

mlogh(ind)=[];

x_alpha(ind)=[];

w(ind)=[];

Elogh=std(logh,[],2);

mU=kb*T*(mlogh-min(mlogh));

sigma2_U=kb*T*Elogh;

% Approximation to a quadratic function, Using linear ftting with weights

c=fit(x_alpha',mlogh,'poly2','weights',w);

%return to original variables after fit

x_alpha=x_alpha*maxbin;

k_pot=2*kb*T*c.p1/maxbin^2;

x_eq=-c.p2/(2*c.p1)*maxbin;

rho0=exp(c.p1/maxbin^2*x_eq^2-c.p3);

cint=confint(c,0.68); %0.68 corresponds to one standard deviation


sigma2_k_pot=2*kb*T/maxbin^2*(cint(2,1)-cint(1,1))/2;



%
disp('...')

disp('Potential analysis using linear fitting')

disp(['k_pot: ' num2str(k_pot*1e6) '+-' num2str(sigma2_k_pot*1e6) ' pN/um']);

