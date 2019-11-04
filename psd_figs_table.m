close all, clear all;

%%  ==============Parameter declaration============


kb=1.38064852e-23; % Boltzmann constant [m^2kg/s^2K]

r=1.03E-6;      % Particle radius [m]



xwi = 400;    % width of the plot square
bx1 = 150;     % extra space at the left
bx2 = 30;     % extra space at the right

Xpix = 3*xwi+3*bx1+3*bx2;  % total

ywi = 300;    % length riquadro con funzione
by1 = 110;     % extra space below
by2 = 50;     % extra space up

Ypix = 1*by1+1*ywi+1*by2;  % larghezza figura in pixel
%number of bins of the histogram, if not set default is 50


%use a subsampled data set
subs=3;

%%  =========Loading selected file============
%load('Data_positions_Fig9_1P6_S.mat')

addpath psd


%Boltzmann constant
kb=1.38064852e-23;




%% plot figures



%creates the figure to do the subplots
figure('Position',[10 20 Xpix Ypix]);
%first figure, probability distribution, exp I

%%
titleI='Experiment I, P=2.3mW';
[k_psd_I, Ek_psd_I, mgamma_psd_I, Egamma_psd_I]=plotsub_psd('Data_positions_Fig9_1P2_S.mat',[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix, titleI, subs);
disp('................')
disp(titleI)
disp('Power spectrum density analysis by linear fitting')
disp(['k_lf: ' num2str(k_psd_I*1e6) '+-' num2str(Ek_psd_I*1e6) 'pN/um'])
disp(['gamma_acf:' num2str(mgamma_psd_I*1e9) '+-'  num2str(Egamma_psd_I*1e9) ' pN ms/um ']);


%%
titleII='Experiment II, P=6.0mW';
[k_psd_II, Ek_psd_II, mgamma_psd_II, Egamma_psd_II]=plotsub_psd('Data_positions_Fig9_1P4_S.mat',[2*bx1+xwi+bx2 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix, titleII, subs);
disp('................')
disp(titleII)
disp('Power spectrum density analysis by linear fitting')
disp(['k_lf: ' num2str(k_psd_II*1e6) '+-' num2str(Ek_psd_II*1e6) 'pN/um'])
disp(['gamma_acf:' num2str(mgamma_psd_II*1e9) '+-'  num2str(Egamma_psd_II*1e9) ' pN ms/um ']);

%%
titleIII='Experiment III, P=9.2mW';
[k_psd_III, Ek_psd_III, mgamma_psd_III, Egamma_psd_III]=plotsub_psd('Data_positions_Fig9_1P6_S.mat',[3*bx1+2*xwi+2*bx2 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix, titleIII, subs);


disp('................')
disp(titleIII)
disp('Power spectrum density analysis by linear fitting')
disp(['k_lf: ' num2str(k_psd_III*1e6) '+-' num2str(Ek_psd_III*1e6) 'pN/um'])
disp(['gamma_acf:' num2str(mgamma_psd_III*1e9) '+-'  num2str(Egamma_psd_III*1e9) ' pN ms/um ']);

