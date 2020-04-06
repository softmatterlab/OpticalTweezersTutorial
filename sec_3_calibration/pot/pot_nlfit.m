function [k_pot, sigma2_k_pot, x_alpha, mrho, sigma2_rho, mU, sigma2_U, rho0, x_eq,  U_0_exp]=pot_nlfit(x,T,varargin)
    %POT_NLFIT   1D implementation of the POTENTIAL METHOD using non-linear fitting.
    %   POT_NLFIT(X,T, P) generates a estimator of the stiffness k_pot
    %   for the potential method using non-linear fitting and fitting
    %   function ft=fittype('sqrt(a/pi)*exp(-a*(x-b)^2)')
    %
    %   INPUTS
    %   x: time series of the position of the particle
    %   T: temperature
    %   P: (optional) number of bins for the histogram, default is 50
    %
    %   OUTPUTS
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

kb = 1.38064852e-23;

%translate everithing to zero
x = x - repmat(mean(x),size(x,1),1);

%default number of bins
P=50;

%user defined number of bins
if nargin>2
    P=varargin{1};
end
[x_alpha, mrho, sigma2_rho, ~, ~, mlogh, sigma2_logh,mU,  sigma2_U, ~, ~, U_0_exp] = prob_dist_energy(x,P, T);
%delete zeros to avoid Inf in weights
sigma2_rho(sigma2_rho==0) = 1;

% weights for fitting
w=1./sigma2_rho.^2;

%in case any other Inf value arises
w(isinf(w)) = 1;

%normalization to avoid "Equation is badly conditioned"
maxbin =max(abs(x_alpha)); 
x_alpha = x_alpha/maxbin; 

%%delete values of zero probability 
ind=find(mrho==Inf);
mrho(ind)=[];
%and all the corresponding values of that data point in 
%%order to keep the same arra size for all variables
x_alpha(ind)=[];
w(ind)=[];
sigma2_rho(ind)=[];
%size(sigma2_rho)
%size(mrho)
mU(ind)=[];
sigma2_U(ind)=[];

% Guess for the initial conditions for the non-linear fitting
a0=(max(mrho)*sqrt(pi))^2;
b0=0;

%Using non-linear fitting with weights
ft=fittype('sqrt(a/pi)*exp(-a*(x-b)^2)');
c=fit(x_alpha',mrho*maxbin,ft,'weights',w/maxbin^2,'StartPoint',[a0*maxbin^2,b0]);

%return to original variables after fit
x_alpha=x_alpha*maxbin;

%stiffness
k_pot=2*kb*T*c.a/maxbin^2; 

%equilibrium position
x_eq=c.b*maxbin;

%0.68 corresponds to one standard deviation
cint=confint(c,0.95);  

%standard deviation squared for the stiffness
sigma2_k_pot=2*kb*T/maxbin^2*(cint(2,1)-cint(1,1))/2;



%estimation for the normalization factor
rho0=sqrt(c.a/pi)*1/maxbin;  

%disp('...')

%disp('Potential analysis using non linear fitting')

%disp(['k_pot: ' num2str(k_pot*1e6) '+-' num2str(sigma2_k_pot*1e6) ' pN/um']);

