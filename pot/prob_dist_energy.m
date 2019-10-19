function [x_alpha, mrho, sigma2_rho, frequencynor, logf, mlogf, sigma_logh,mU,  sigma2_U, mloghist, Eloghist]=prob_dist_energy(x,P, T)
% PROB_DIST(x, T, P) defines the probability distribution for a set of x
% positions and P number of bins
%
%
%
% INPUTS
% x     :   positions
% P     :   number of bins
%
%OUTPUTS
%x_alpha: central position of the bins
%mrho: probability distribution
%sigma2_rho: standard deviation of the probability distribution
%frequency: counts for each experiment


kb=1.38064852e-23;

%obtain number of experiments
[~,Nexp]=size(x);
%defines the bins edges
edges=linspace(min(x(:)),max(x(:)),P);

%defines the lenght of bins
dx=edges(2)-edges(1);

%defines central position of the bins
x_alpha=(edges(2:end)+edges(1:end-1))/2;

%define the histogram
for j=1:Nexp
    
    xx=x(:,j);

    [frequency(:,j)]=histcounts(xx,edges);
    
end

% after normalization it becomes the probability density distribution
frequencynor=frequency./(sum(frequency,1)*dx);
disp(size(frequencynor))
%mean probability  density distribution
mrho=mean(frequencynor,2);
disp(size(mrho));
%logaritm of the histogram 
loghist=-log(frequency);
mloghist=mean(loghist, 2);
Eloghist=std(loghist,[], 2);

%standard deviation squared of probability density distribution
sigma2_rho=std(frequencynor,[],2);
disp(size(sigma2_rho))

%log of the frequency
logf=-log(frequencynor);

%mean log of the frequency
mlogf=mean(logf,2);

%standard deviation of the log of the frequency
sigma_logh=std(logf,[],2);

%mean value of the potential energy
mU=kb*T*(mlogf-min(mlogf));

%standard deviation of the potential energy
sigma2_U=kb*T*sigma_logh;

end