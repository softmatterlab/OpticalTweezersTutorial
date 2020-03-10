% PSD analysis for Nexp

% Initialization of the workspace
clear;

close all;

addpath forma



load('Data_positions_Fig9_1P2_S.mat')


x = x - repmat(mean(x),size(x,1),1);

gamma=6*pi*eta*a;

kb=1.38064852e-23;

D0=kb*T/gamma;

subs=10; %use a subsampled data set

[N,Nexp]=size(x);

[fc_forma,D_forma,Efcsol,EDsol]=forma1d(x(1:subs:N,:),dt*subs);

mD_forma=mean(D_forma);

ED_forma=std(D_forma);

gamma_forma=kb*T./D_forma;

mgamma_forma=mean(gamma_forma);

Egamma_forma=std(gamma_forma);

k_forma=2*pi*gamma*fc_forma;

mk_forma=mean(k_forma);

Ek_forma=std(k_forma);

% estimation of k using the estimated gamma
k2_forma=2*pi*gamma_forma.*fc_forma;

mk2_forma=mean(k2_forma);

Ek2_forma=std(k2_forma);

Eksol=2*pi*gamma_forma.*Efcsol;



disp('...')

disp('FORMA analysis')

disp(['D: ' num2str(mD_forma) '+-' num2str(ED_forma) ' m^2/s'])

disp(['D0: ' num2str(D0) ' m^2/s']);

disp(['gamma: ' num2str(mgamma_forma) '+-' num2str(Egamma_forma) ' Ns/m'])

disp(['gamma0: ' num2str(gamma) ' Ns/m'])

disp(['k: ' num2str(mk_forma) '+-' num2str(Ek_forma) ' N/m'])

disp(['k2: ' num2str(mk2_forma) '+-' num2str(Ek2_forma) ' N/m'])

disp(['fc: ' num2str(mean(fc_forma)) ' Hz'])

disp(['Esol:' num2str(Eksol)])

