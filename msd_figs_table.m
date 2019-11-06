close all;
clear all;

%%  ==============Parameter declaration============


kB=1.38064852e-23; % Boltzmann constant [m^2kg/s^2K]



xwi = 400;    % width of the plot square
bx1 = 120;     % extra space at the left
bx2 = 20;     % extra space at the right

Xpix = 3*xwi+bx1+3*bx2;  % total

ywi = 300;    % length riquadro con funzione
by1 = 110;     % extra space below
by2 = 70;     % extra space up

Ypix = 1*by1+1*ywi+1*by2;  % larghezza figura in pixel
T=293.15;

%%  =========Loading selected file============
%load('Data_positions_Fig9_1P6_S.mat')

addpath msd
addpath wlsice


%Boltzmann constant
kB=1.38064852e-23;




%% plot figures



%creates the figure to do the subplots
figure('Position',[10 20 Xpix Ypix]);
%first figure, probability distribution, exp I

%%
subs=1;
 maxlag=1000;

titleI='Experiment I, P=2.3mW';
[k_msd_I,sigma_k_msd_I,  gamma_msd_I, sigma_gamma_msd_I,D_msd_I,sigma_D_msd_I, tau0_I]=plotsub_msd('Data_positions_Fig9_1P2_S.mat',[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix,  titleI, T, 20, 50,0.8,3,5);

                        %                                            plotsub_msd(filename,                           positioninthefig1,                      title1, T, subs, maxlag,partau0,ytau)
disp('................')
disp(titleI)
disp('Autocorrelation function analysis by linear fitting')
 
disp(['k_msd: ' num2str(k_msd_I*1e6) '+-' num2str(sigma_k_msd_I*1e6) ' pN/um'])
 

disp(['D_msd: ' num2str(D_msd_I*1e12) '+-' num2str(sigma_D_msd_I*1e12) ' um^2/s'])
 
disp(['gamma_acf:' num2str(gamma_msd_I*1e9) '+-'  num2str(sigma_gamma_msd_I*1e9) ' pN s/um ']);

disp(['tau_0:' num2str(tau0_I*1e3) ' ms']);
disp('................')
%%
titleII='Experiment II, P=6.0mW';
[k_msd_II,sigma2_k_msd_II,  gamma_msd_II, sigma2_gamma_msd_II,D_msd_II, sigma_D_msd_II, tau0_II]=plotsub_msd('Data_positions_Fig9_1P4_S.mat',[bx1+bx2+xwi 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix, titleII, T, 17, 70,1.15,3,1);

disp('................')
disp(titleII)
disp('Autocorrelation function analysis by linear fitting')
 
disp(['k_msd: ' num2str(k_msd_II*1e6) '+-' num2str(sigma2_k_msd_II*1e6) ' pN/um'])
 
disp(['D_msd: ' num2str(D_msd_II*1e12) '+-' num2str(sigma_D_msd_II*1e12) ' um^2/s']) 
disp(['gamma_acf:' num2str(gamma_msd_II*1e9) '+-'  num2str(sigma2_gamma_msd_II*1e9) ' pN s/um ']);

disp(['tau_0:' num2str(tau0_II*1e3) ' ms']);
disp('................')
%%
titleIII='Experiment III, P=9.2mW';
[k_msd_III,sigma2_k_msd_III,  gamma_msd_III, sigma2_gamma_msd_III, D_msd_III, sigma_D_msd_III, tau0_III]=plotsub_msd('Data_positions_Fig9_1P6_S.mat',[bx1+2*bx2+2*xwi 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix, titleIII, T, 10, 100,1.2,3,1);


disp('................')
disp(titleIII)
disp('Autocorrelation function analysis by linear fitting')
 
disp(['k_msd: ' num2str(k_msd_III*1e6) '+-' num2str(sigma2_k_msd_III*1e6) ' pN/um'])
 
disp(['D_msd: ' num2str(D_msd_III*1e12) '+-' num2str(sigma_D_msd_III*1e12) ' um^2/s']) 
disp(['gamma_acf:' num2str(gamma_msd_III*1e9) '+-'  num2str(sigma2_gamma_msd_III*1e9) ' pN s/um ']);

disp(['tau_0:' num2str(tau0_III*1e3) ' ms']);
disp('................')
