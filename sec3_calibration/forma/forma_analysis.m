% PSD analysis for Nexp

% Initialization of the workspace
clear;

close all;

addpath ../data/
addpath ../statistics_func/

%%
%first experiment
disp('first experiment')
load('Data_x_positions_Exp_I.mat')


x = x - repmat(mean(x),size(x,1),1);
%in this case the whole time series is used as one experiment but it can be
%done with many experiments
xl=reshape(x, [size(x,1)*size(x,2),1 ]);

gamma=6*pi*eta*a;

kb=1.38064852e-23;

D0=kb*T/gamma;

subs=1; %use a subsampled data set
nsubs=3;
[N,Nexp]=size(x);



[fc_forma,D_forma,sigma_fcforma,sigma_D_forma]=forma1d(x(1:subs:end, 1),dt*subs, nsubs);

gamma_forma=kb*T./D_forma;
sigma_gamma_forma=kb*T./D_forma.^2.*sigma_D_forma;


k_forma=gamma_forma.*fc_forma;
sigma_k_forma=gamma_forma.*sigma_fcforma+fc_forma.*sigma_gamma_forma;




disp('................')

disp('FORMA analysis for a single experiment using formula 41')
disp(['k_forma: ' num2str(k_forma*1e6) '+-' num2str(sigma_k_forma*1e6) 'p Nu/m'])

disp(['D_forma: ' num2str(D_forma*1e12) '+-' num2str(sigma_D_forma*1e12) 'u m^2/s'])

disp(['gamma_forma: ' num2str(1e9*gamma_forma) '+-' num2str(1e9*sigma_gamma_forma) ' pNms/um'])


disp('................')


%[fc_forma,D_forma,sigma_fcforma,sigma_D_forma]=forma1d(xl(1:subs:end),dt*subs, nsubs);
for jj=1:Nexp
[fc_forma(jj),D_sforma(jj),sigma_fcforma,sigma_D_forma]=forma1d(x(:, jj),dt*subs, nsubs);

gamma_sforma(jj)=kb*T./D_sforma(jj);
%sigma_gamma_forma=kb*T./D_sforma.^2.*sigma_D_forma;


k_forma_s(jj)=gamma_sforma(jj).*fc_forma(jj);
%sigma_k_forma_s(jj)=gamma_forma.*sigma_fcforma+fc_forma.*sigma_gamma_forma;


end 

disp('................')

disp('FORMA analysis averaging all the experiments')


disp(['k_forma: ' num2str(mean(k_forma_s)*1e6) '+-' num2str(std(k_forma_s)*1e6) 'p Nu/m'])

 disp(['D_forma: ' num2str(mean(D_sforma)*1e12) '+-' num2str(std(D_sforma)*1e12) 'u m^2/s'])
 
 disp(['gamma_forma: ' num2str(1e9*mean(gamma_sforma)) '+-' num2str(1e9*std(gamma_sforma)) ' pNms/um'])
% 
% 
disp('................')


k=mean(k_forma_s)*1e6;
dk=std(k_forma_s)*1e6;
D=mean(D_sforma)*1e12;
dD=std(D_sforma)*1e12;
g=1e9*mean(gamma_sforma);
dg=1e9*std(gamma_sforma);

[v1, dv1, sig]=round_significance(k, dk);

[v1, dv1, sig]=round_significance(D, dD);

[v1, dv1, sig]=round_significance(g, dg);



[fc_forma_l,D_forma_l,sigma_fcforma_l,sigma_D_forma_l]=forma1d(xl(1:subs:end),dt*subs, nsubs);



gamma_forma_l=kb*T./D_forma_l;
sigma_gamma_forma_l=kb*T./D_forma_l.^2.*sigma_D_forma_l;


k_forma_l=gamma_forma_l.*fc_forma_l;
sigma_k_forma_l=gamma_forma_l.*sigma_fcforma_l+fc_forma_l.*sigma_gamma_forma_l;




disp('................')

disp('FORMA analysis for all the experiments as a single experiment and using f.41')
disp(['k_forma: ' num2str(k_forma_l*1e6) '+-' num2str(sigma_k_forma_l*1e6) 'p Nu/m'])

disp(['D_forma: ' num2str(D_forma_l*1e12) '+-' num2str(sigma_D_forma_l*1e12) 'u m^2/s'])

disp(['gamma_forma: ' num2str(1e9*gamma_forma_l) '+-' num2str(1e9*sigma_gamma_forma_l) ' pNms/um'])


disp('................')


%%
%second experiment
disp('second experiment')


load('Data_x_positions_Exp_II.mat')


x = x - repmat(mean(x),size(x,1),1);
%in this case the whole time series is used as one experiment but it can be
%done with many experiments
xl=reshape(x, [size(x,1)*size(x,2),1 ]);

gamma=6*pi*eta*a;

kb=1.38064852e-23;

D0=kb*T/gamma;

subs=1; %use a subsampled data set
nsubs=3;
[N,Nexp]=size(x);

%[fc_forma,D_forma,sigma_fcforma,sigma_D_forma]=forma1d(xl(1:subs:end),dt*subs, nsubs);

[fc_forma,D_forma,sigma_fcforma,sigma_D_forma]=forma1d(x(1:subs:end, 1),dt*subs, nsubs);

gamma_forma=kb*T./D_forma;
sigma_gamma_forma=kb*T./D_forma.^2.*sigma_D_forma;


k_forma=gamma_forma.*fc_forma;
sigma_k_forma=gamma_forma.*sigma_fcforma+fc_forma.*sigma_gamma_forma;




disp('................')

disp('FORMA analysis for a single experiment using formula 41')
disp(['k_forma: ' num2str(k_forma*1e6) '+-' num2str(sigma_k_forma*1e6) 'p Nu/m'])

disp(['D_forma: ' num2str(D_forma*1e12) '+-' num2str(sigma_D_forma*1e12) 'u m^2/s'])

disp(['gamma_forma: ' num2str(1e9*gamma_forma) '+-' num2str(1e9*sigma_gamma_forma) ' pNms/um'])


disp('................')


%[fc_forma,D_forma,sigma_fcforma,sigma_D_forma]=forma1d(xl(1:subs:end),dt*subs, nsubs);
for jj=1:Nexp
[fc_forma(jj),D_sforma(jj),sigma_fcforma,sigma_D_forma]=forma1d(x(:, jj),dt*subs, nsubs);

gamma_sforma(jj)=kb*T./D_sforma(jj);
%sigma_gamma_forma=kb*T./D_sforma.^2.*sigma_D_forma;


k_forma_s(jj)=gamma_sforma(jj).*fc_forma(jj);
%sigma_k_forma_s(jj)=gamma_forma.*sigma_fcforma+fc_forma.*sigma_gamma_forma;


end 

disp('................')

disp('FORMA analysis averaging all the experiments')


disp(['k_forma: ' num2str(mean(k_forma_s)*1e6) '+-' num2str(std(k_forma_s)*1e6) 'p Nu/m'])

 disp(['D_forma: ' num2str(mean(D_sforma)*1e12) '+-' num2str(std(D_sforma)*1e12) 'u m^2/s'])
 
 disp(['gamma_forma: ' num2str(1e9*mean(gamma_sforma)) '+-' num2str(1e9*std(gamma_sforma)) ' pNms/um'])
% 
% 
disp('................')


k=mean(k_forma_s)*1e6;
dk=std(k_forma_s)*1e6;
D=mean(D_sforma)*1e12;
dD=std(D_sforma)*1e12;
g=1e9*mean(gamma_sforma);
dg=1e9*std(gamma_sforma);

[v1, dv1, sig]=round_significance(k, dk);

[v1, dv1, sig]=round_significance(D, dD);

[v1, dv1, sig]=round_significance(g, dg);



[fc_forma_l,D_forma_l,sigma_fcforma_l,sigma_D_forma_l]=forma1d(xl(1:subs:end),dt*subs, nsubs);



gamma_forma_l=kb*T./D_forma_l;
sigma_gamma_forma_l=kb*T./D_forma_l.^2.*sigma_D_forma_l;


k_forma_l=gamma_forma_l.*fc_forma_l;
sigma_k_forma_l=gamma_forma_l.*sigma_fcforma_l+fc_forma_l.*sigma_gamma_forma_l;




disp('................')

disp('FORMA analysis for all the experiments as a single experiment and using f.41')
disp(['k_forma: ' num2str(k_forma_l*1e6) '+-' num2str(sigma_k_forma_l*1e6) 'p Nu/m'])

disp(['D_forma: ' num2str(D_forma_l*1e12) '+-' num2str(sigma_D_forma_l*1e12) 'u m^2/s'])

disp(['gamma_forma: ' num2str(1e9*gamma_forma_l) '+-' num2str(1e9*sigma_gamma_forma_l) ' pNms/um'])


disp('................')




%%
%third experiment
disp('third experiment')

load('Data_x_positions_Exp_III.mat')


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

%[fc_forma,D_forma,sigma_fcforma,sigma_D_forma]=forma1d(xl(1:subs:end),dt*subs, nsubs);

[fc_forma,D_forma,sigma_fcforma,sigma_D_forma]=forma1d(x(1:subs:end, 1),dt*subs, nsubs);

gamma_forma=kb*T./D_forma;
sigma_gamma_forma=kb*T./D_forma.^2.*sigma_D_forma;


k_forma=gamma_forma.*fc_forma;
sigma_k_forma=gamma_forma.*sigma_fcforma+fc_forma.*sigma_gamma_forma;




disp('................')

disp('FORMA analysis for a single experiment using formula 41')
disp(['k_forma: ' num2str(k_forma*1e6) '+-' num2str(sigma_k_forma*1e6) 'p Nu/m'])

disp(['D_forma: ' num2str(D_forma*1e12) '+-' num2str(sigma_D_forma*1e12) 'u m^2/s'])

disp(['gamma_forma: ' num2str(1e9*gamma_forma) '+-' num2str(1e9*sigma_gamma_forma) ' pNms/um'])


disp('................')


%[fc_forma,D_forma,sigma_fcforma,sigma_D_forma]=forma1d(xl(1:subs:end),dt*subs, nsubs);
for jj=1:Nexp
[fc_forma(jj),D_sforma(jj),sigma_fcforma,sigma_D_forma]=forma1d(x(:, jj),dt*subs, nsubs);

gamma_sforma(jj)=kb*T./D_sforma(jj);
%sigma_gamma_forma=kb*T./D_sforma.^2.*sigma_D_forma;


k_forma_s(jj)=gamma_sforma(jj).*fc_forma(jj);
%sigma_k_forma_s(jj)=gamma_forma.*sigma_fcforma+fc_forma.*sigma_gamma_forma;


end 

disp('................')

disp('FORMA analysis averaging all the experiments')


disp(['k_forma: ' num2str(mean(k_forma_s)*1e6) '+-' num2str(std(k_forma_s)*1e6) 'p Nu/m'])

 disp(['D_forma: ' num2str(mean(D_sforma)*1e12) '+-' num2str(std(D_sforma)*1e12) 'u m^2/s'])
 
 disp(['gamma_forma: ' num2str(1e9*mean(gamma_sforma)) '+-' num2str(1e9*std(gamma_sforma)) ' pNms/um'])
% 
% 
disp('................')


k=mean(k_forma_s)*1e6;
dk=std(k_forma_s)*1e6;
D=mean(D_sforma)*1e12;
dD=std(D_sforma)*1e12;
g=1e9*mean(gamma_sforma);
dg=1e9*std(gamma_sforma);

[v1, dv1, sig]=round_significance(k, dk);


[v1, dv1, sig]=round_significance(D, dD);


[v1, dv1, sig]=round_significance(g, dg);




[fc_forma_l,D_forma_l,sigma_fcforma_l,sigma_D_forma_l]=forma1d(xl(1:subs:end),dt*subs, nsubs);



gamma_forma_l=kb*T./D_forma_l;
sigma_gamma_forma_l=kb*T./D_forma_l.^2.*sigma_D_forma_l;


k_forma_l=gamma_forma_l.*fc_forma_l;
sigma_k_forma_l=gamma_forma_l.*sigma_fcforma_l+fc_forma_l.*sigma_gamma_forma_l;




disp('................')

disp('FORMA analysis for all the experiments as a single experiment and using f.41')
disp(['k_forma: ' num2str(k_forma_l*1e6) '+-' num2str(sigma_k_forma_l*1e6) 'p Nu/m'])

disp(['D_forma: ' num2str(D_forma_l*1e12) '+-' num2str(sigma_D_forma_l*1e12) 'u m^2/s'])

disp(['gamma_forma: ' num2str(1e9*gamma_forma_l) '+-' num2str(1e9*sigma_gamma_forma_l) ' pNms/um'])


disp('................')

