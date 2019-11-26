close all, clear all;

%%  ==============Parameter declaration============


kB=1.38e-23; % Boltzmann constant [m^2kg/s^2K]
T=300;  % Temperature [K]
r=1.03E-6;      % Particle radius [m]
v=0.00002414*10^(247.8/(-140+T));  % Water viscosity [Pa*s]
gamma=pi*6*r*v; %[m*Pa*s]

xwi = 400;    % width of the plot square
bx1 = 120;     % extra space at the left
bx2 = 20;     % extra space at the right

Xpix = 3*xwi+bx1+3*bx2;  % total

ywi = 300;    % length riquadro con funzione
by1 = 110;     % extra space below
by2 = 70;     % extra space up

Ypix = 1*by1+1*ywi+1*by2;  % larghezza figura in pixel
%number of bins of the histogram, if not set default is 50
P=50; 

%use a subsampled data set
subs=1;

%%  =========Loading selected file============
%load('Data_positions_Fig9_1P6_S.mat')

addpath acf


%Boltzmann constant
kb=1.38064852e-23;




%% plot figures



%creates the figure to do the subplots
figure('Position',[10 20 Xpix Ypix]);
%first figure, probability distribution, exp I

%%
titleI='Experiment I, P=2.3mW';
[k_acf_lf_I,sigma_k_acf_lf_I,D_acf_lf_I,sigma_D_acf_lf_I,gamma_acf_lf_I, sigma_gamma_acf_lf_I , k_acf_nl_I, sigma_k_acf_nl_I,...
D_acf_nl_I, sigma_D_acf_nl_I,gamma_acf_nl_I, sigma_gamma_acf_nl_I, tau0_exp_lf_I,sigma_tau0_exp_lf_I,  tau0_exp_nl_I, sigma_tau0_exp_nl_I, Xpix, Ypix]=...
plotsub_acf('Data_positions_Fig9_1P2_S.mat',[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix, titleI, T, subs,1.1,3.21,5);
 disp('................')
 disp(titleI)
 disp('Autocorrelation function analysis by linear fitting')
 
disp(['k_acf: ' num2str(k_acf_lf_I*1e6) '+-' num2str(sigma_k_acf_lf_I*1e6) ' pN/um'])
 
disp(['D_acf: ' num2str(D_acf_lf_I*1e12) '+-' num2str(sigma_D_acf_lf_I*1e12) ' um^2/s'])
 
disp(['gamma_acf:' num2str(gamma_acf_lf_I*1e9) '+-'  num2str(sigma_gamma_acf_lf_I*1e9) ' pN s/um ']);
disp(['tau_0:' num2str(tau0_exp_lf_I*1e3)  '+-' num2str(sigma_tau0_exp_lf_I*1e3)  ' ms']);
 
disp('................')

disp('Autocorrelation function analysis by non-linear fitting')
 
disp(['k_acf: ' num2str(k_acf_nl_I*1e6) '+-' num2str(sigma_k_acf_nl_I*1e6) ' pN/um'])
 
disp(['D_acf: ' num2str(D_acf_nl_I*1e12) '+-' num2str(sigma_D_acf_nl_I*1e12) ' um^2/s'])
 
disp(['gamma_acf:' num2str(gamma_acf_nl_I*1e9) '+-'  num2str(sigma_gamma_acf_nl_I*1e9) ' pN ms/um ']);
disp(['tau_0:' num2str(tau0_exp_nl_I*1e3) '+-' num2str(sigma_tau0_exp_nl_I*1e3)   ' ms']);
 
disp('................')
%%
titleII='Experiment II, P=6.0mW';
[k_acf_lf_II,sigma_k_acf_lf_II,D_acf_lf_II,sigma_D_acf_lf_II,gamma_acf_lf_II, sigma_gamma_acf_lf_II , k_acf_nl_II, sigma_k_acf_nl_II,...
 D_acf_nl_II, sigma_D_acf_nl_II,gamma_acf_nl_II, sigma_gamma_acf_nl_II, tau0_exp_lf_II, sigma_tau0_exp_lf_II , tau0_exp_nl_II,sigma_tau0_exp_nl_II ]=....
plotsub_acf('Data_positions_Fig9_1P4_S.mat',[bx1+bx2+xwi 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix, titleII, T, subs,1.1,8.2,0);

disp('................')
disp(titleII)
disp('Autocorrelation function analysis by linear fitting')
 
disp(['k_acf: ' num2str(k_acf_lf_II*1e6) '+-' num2str(sigma_k_acf_lf_II*1e6) ' pN/um'])
 
disp(['D_acf: ' num2str(D_acf_lf_II*1e12) '+-' num2str(sigma_D_acf_lf_II*1e12) ' um^2/s'])
 
disp(['gamma_acf:' num2str(gamma_acf_lf_II*1e9) '+-'  num2str(sigma_gamma_acf_lf_II*1e9) ' pN ms/um ']);
disp(['tau_0:' num2str(tau0_exp_lf_II*1e3)  '+-' num2str(sigma_tau0_exp_lf_II*1e3)  ' ms']);
 
disp('................')

disp('Autocorrelation function analysis by non-linear fitting')
 
disp(['k_acf: ' num2str(k_acf_nl_II*1e6) '+-' num2str(sigma_k_acf_nl_II*1e6) ' pN/um'])
 
disp(['D_acf: ' num2str(D_acf_nl_II*1e12) '+-' num2str(sigma_D_acf_nl_II*1e12) ' um^2/s'])
 
disp(['gamma_acf:' num2str(gamma_acf_nl_II*1e9) '+-'  num2str(sigma_gamma_acf_nl_II*1e9) ' pN ms/um ']);
disp(['tau_0:' num2str(tau0_exp_nl_II*1e3)  '+-' num2str(sigma_tau0_exp_nl_II*1e3)  ' ms']);
 
disp('................')
%%
titleIII='Experiment III, P=9.2mW';
[k_acf_lf_III,sigma_k_acf_lf_III,D_acf_lf_III,sigma_D_acf_lf_III,gamma_acf_lf_III, sigma_gamma_acf_lf_III , k_acf_nl_III, sigma_k_acf_nl_III,...
 D_acf_nl_III, sigma_D_acf_nl_III,gamma_acf_nl_III, sigma_gamma_acf_nl_III, tau0_exp_lf_III, sigma_tau0_exp_lf_III , tau0_exp_nl_III,sigma_tau0_exp_nl_III]=...
plotsub_acf('Data_positions_Fig9_1P6_S.mat',[bx1+2*bx2+2*xwi 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix, titleIII, T, subs, 1.1,13,0);

disp('................')
disp(titleIII)
disp('Autocorrelation function analysis by linear fitting')
 
disp(['k_acf: ' num2str(k_acf_lf_III*1e6) '+-' num2str(sigma_k_acf_lf_III*1e6) ' pN/um'])
 
disp(['D_acf: ' num2str(D_acf_lf_III*1e12) '+-' num2str(sigma_D_acf_lf_III*1e12) ' um^2/s'])
 
disp(['gamma_acf:' num2str(gamma_acf_lf_III*1e9) '+-'  num2str(sigma_gamma_acf_lf_III*1e9) ' pN s/um ']);
disp(['tau_0:' num2str(tau0_exp_lf_III*1e3)  '+-' num2str(sigma_tau0_exp_lf_III*1e3)  ' ms']);
 
disp('................')

disp('Autocorrelation function analysis by non-linear fitting')
 
disp(['k_acf: ' num2str(k_acf_nl_III*1e6) '+-' num2str(sigma_k_acf_nl_III*1e6) ' pN/um'])
 
disp(['D_acf: ' num2str(D_acf_nl_III*1e12) '+-' num2str(sigma_D_acf_nl_III*1e12) ' um^2/s'])
 
disp(['gamma_acf:' num2str(gamma_acf_nl_III*1e9) '+-'  num2str(sigma_gamma_acf_nl_III*1e9) ' pN ms/um ']);
disp(['tau_0:' num2str(tau0_exp_nl_III*1e3)  '+-' num2str(sigma_tau0_exp_nl_III*1e3)  ' ms']);
 
disp('................')
