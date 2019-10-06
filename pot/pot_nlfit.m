function [k_pot, sigma2_k_pot, x_alpha, mrho, sigma2_rho, mU, sigma2_U, rho0, x_eq]=pot_nlfit(x,T,varargin)
    %POT_NLFIT   1D implementation of the POTENTIAL METHOD using non-linear fitting.
    %   POT_NLFIT(X,T, nb) generates a estimator of the stiffness k_pot
    %   for the potential method using non-linear fitting. 
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

ind=find(mrho==Inf);

mrho(ind)=[];

x_alpha(ind)=[];


% Guess for the initial conditions for the non-linear fitting

a0=(max(mrho)*sqrt(pi))^2;

b0=0;




%Using non-linear ftting with weights

ft=fittype('sqrt(a/pi)*exp(-a*(x-b)^2)');

c=fit(x_alpha',mrho*maxbin,ft,'weights',w/maxbin^2,'StartPoint',[a0*maxbin^2,b0]);

%return to original variables after fit

x_alpha=x_alpha*maxbin;

k_pot=2*kb*T*c.a/maxbin^2; %stiffness

x_eq=c.b*maxbin^2;

cint=confint(c,0.68);  %0.68 corresponds to one standard deviation

sigma2_k_pot=2*kb*T/maxbin^2*(cint(2,1)-cint(1,1))/2;

%estimate energy potential
logh=log(frequency);

mlogh=mean(logh,2);

Delta_logh=std(logh,[],2);

mU=-kb*T*(mlogh-max(mlogh));

sigma2_U=-kb*T*Delta_logh;

rho0=sqrt(c.a/pi)*1/maxbin;  %estimation for the normalization factor


disp('...')

disp('Potential analysis using non linear fitting')

disp(['k_pot: ' num2str(k_pot*1e6) '+-' num2str(sigma2_k_pot*1e6) ' pN/um']);

