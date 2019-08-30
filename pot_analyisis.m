% Potential analysis for Nexp

% Initialization of the workspace
clear;

close all;

addpath pot



% Load data file

load(['Fig9_1P4_S.mat']);

Vx=Vx*Sdx;

kb=1.38064852e-23;

subs=1; %use a subsampled data set


[xbins,mU,EU,k_pot,Ek_pot,mhist,Ehist,h0,x0]=pot_nlfit(Vx(1:subs:end,:),T,50);

%[xbins,mU,EU,k_pot,Ek_pot,mhist,Ehist,h0,x0]=pot_nlfit(Vx(1:subs:end,:),T,50);


disp(['tau_0: ' num2str(6*pi*eta*a/k_pot)])
disp(['dt: ' num2str(dt*subs)])

%save([nfile filesep datafile '_fit_pot_subs' num2str(subs)  '.mat'],'xbins','mU','EU','k_pot','Ek_pot','mhist','Ehist','h0','x0');
 
%savePDF(gcf,[nfile filesep datafile '_pot_subs' num2str(subs) '.pdf'])

