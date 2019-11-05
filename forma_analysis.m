% PSD analysis for Nexp

% Initialization of the workspace
clear;

close all;

addpath forma



load('Data_positions_Fig9_1P6_S.mat')


x = x - repmat(mean(x),size(x,1),1);
%in this case the whole time series is used as one experiment but it can be
%done wioth many experiments
xl=reshape(x, [size(x,1)*size(x,2),1 ]);

gamma=6*pi*eta*a;

kb=1.38064852e-23;

D0=kb*T/gamma;

subs=1; %use a subsampled data set
nsubs=3;
[N,Nexp]=size(x);

[fc_forma,D_forma,sigma_fcforma,sigma_D_forma]=forma1d(xl(1:subs:end),dt*subs, nsubs);


gamma_forma=kb*T./D_forma;
sigma_gamma_forma=kb*T./D_forma^2*sigma_D_forma;


k_forma=gamma_forma*fc_forma;
sigma_k_forma=gamma_forma.*sigma_fcforma+fc_forma*sigma_gamma_forma;




disp('................')

disp('FORMA analysis')
disp(['k_forma: ' num2str(k_forma*1e6) '+-' num2str(sigma_k_forma*1e6) 'p Nu/m'])

disp(['D_forma: ' num2str(D_forma*1e12) '+-' num2str(sigma_D_forma*1e12) 'u m^2/s'])

disp(['gamma_forma: ' num2str(1e9*gamma_forma) '+-' num2str(1e9*sigma_gamma_forma) ' pNms/um'])


disp('................')


