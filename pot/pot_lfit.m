function [k_pot, sigma2_k_pot, x_alpha, mrho, sigma2_rho, mU, sigma2_U, rho0, x_eq, U_0_exp]=pot_lfit(x,T,varargin)
    %POT_LFIT   1D implementation of the POTENTIAL METHOD using linear fitting.
    %   POT_LFIT(X,T, P) generates a estimator of the stiffness k_pot
    %   for the potential method using linear fitting
    %   c=fit(x_alpha',mlogh,'poly2','weights',w)  corresponding to 
    %   c(x) = p1*x^2 + p2*x + p3
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
 
%translate everithing to zero
x = x - repmat(mean(x),size(x,1),1);
 
%default number of bins
P=50;
 
%user defined number of bins
if nargin>2
    P=varargin{1};
end
 
[x_alpha, mrho, sigma2_rho, ~, ~, mlogf, sigma_logh,mU,  sigma2_U, ~, ~, U_0_exp]=prob_dist_energy(x,P, T);
 
%delete zeros to avoid Inf in weights
sigma2_rho(sigma2_rho==0)=1;
 
 
 
 
%delete values with Inf in mlogh
ind=find(isinf(mlogf));
mlogf(ind)=[];
%and all the corresponding values of that data point in 
%%order to keep the same arra size for all variables
x_alpha(ind)=[];
 
sigma2_rho(ind)=[];
mU(ind)=[];
sigma2_U(ind)=[];
mrho(ind)=[];
 
 
% weights for fitting
w=1./sigma_logh.^2;
 
%in case any other Inf value arises
w(isinf(w))=1;
try w(ind)=[];
    w(ind)=[];
catch 
    warning('The value of the weight was already erased');
    
end
 
%normalization to avoid "Equation is badly conditioned"
maxbin=x_alpha(end);
x_alpha=x_alpha/maxbin; 
 
%Approximation to a quadratic function, Using linear fitting with weights
c=fit(x_alpha',mlogf,'poly2','weights',w);
 
 
%return to original variables after fit
x_alpha=x_alpha*maxbin;
 
%stiffness
k_pot=2*kb*T*c.p1/maxbin^2;
 
%equilibrium position
x_eq=-c.p2/(2*c.p1)*maxbin;
 
%estimation for the normalization factor
rho0=exp(c.p1/maxbin^2*x_eq^2-c.p3);
 
%0.68 corresponds to one standard deviation
cint=confint(c,0.95); 
U_0_exp=min(mlogf); %in kbT units
%standard deviation squared for the stiffness
sigma2_k_pot=2*kb*T/maxbin^2*(cint(2,1)-cint(1,1))/2;
 
 
 
%
%disp('...')
 
%disp('Potential analysis using linear fitting')
 
%disp(['k_pot: ' num2str(k_pot*1e6) '+-' num2str(sigma2_k_pot*1e6) ' pN/um']);
 

