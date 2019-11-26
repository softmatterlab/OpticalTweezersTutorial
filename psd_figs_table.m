close all, clear all;

%%  ==============Parameter declaration============


kB=1.38064852e-23; % Boltzmann constant [m^2kg/s^2K]

r=1.03E-6;      % Particle radius [m]



xwi = 400;    % width of the plot square
bx1 = 140;     % extra space at the left
bx2 = 20;     % extra space at the right

Xpix = 3*xwi+bx1+3*bx2;  % total

ywi = 300;    % length riquadro con funzione
by1 = 110;     % extra space below
by2 = 70;     % extra space up

Ypix = 1*by1+1*ywi+1*by2;  % larghezza figura in pixel
%number of bins of the histogram, if not set default is 50


%use a subsampled data set
subs=1;

%%  =========Loading selected file============
%load('Data_positions_Fig9_1P6_S.mat')

addpath psd





%% plot figures

T=293.15

%creates the figure to do the subplots
figure('Position',[10 20 Xpix Ypix]);
%first figure, probability distribution, exp I

%%
titleI='Experiment I, P=2.3mW';
[k_psd_I, sigma_k_psd_I, gamma_psd_I, sigma_gamma_psd_I, D_psd_I, sigma_D_psd_I]=plotsub_psd('Data_positions_Fig9_1P2_S.mat',[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix, titleI, T, subs,5);
disp('................')
disp(titleI)
disp('Power spectrum density analysis by linear fitting')

disp(['k_psd: ' num2str(k_psd_I*1e6) '+-' num2str(sigma_k_psd_I*1e6) 'pN/um'])
[k1, dk1, sig]=round_significance(k_psd_I*1e6, sigma_k_psd_I*1e6);
disp(['k_psd: ' k1 '\pm' dk1 ])

disp(['gamma_psd:' num2str(gamma_psd_I*1e9) '+-'  num2str(sigma_gamma_psd_I*1e9) ' pN ms/um ']);
[k1, dk1, sig]=round_significance(gamma_psd_I*1e9, sigma_gamma_psd_I*1e9);
disp(['gamma_psd:' k1 '\pm' dk1])

disp(['D_psd: ' num2str(D_psd_I*1e12) '+-' num2str(sigma_D_psd_I*1e12) ' um^2/s'])
[k1, dk1, sig]=round_significance(D_psd_I*1e12, sigma_D_psd_I*1e12);
disp(['D_psd: ' k1 '\pm' dk1])

disp('................' )

%%
titleII='Experiment II, P=6.0mW';
[k_psd_II, sigma_k_psd_II, gamma_psd_II, sigma_gamma_psd_II, D_psd_II, sigma_D_psd_II]=plotsub_psd('Data_positions_Fig9_1P4_S.mat',[bx1+bx2+xwi 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix, titleII, T, subs,1);
disp('................')
disp(titleII)
disp('Power spectrum density analysis by linear fitting')
disp(['k_psd: ' num2str(k_psd_II*1e6) '+-' num2str(sigma_k_psd_II*1e6) 'pN/um'])
[k1, dk1, sig]=round_significance(k_psd_II*1e6, sigma_k_psd_II*1e6);
disp(['k_psd: ' k1 '\pm' dk1 ])

disp(['gamma_psd:' num2str(gamma_psd_II*1e9) '+-'  num2str(sigma_gamma_psd_II*1e9) ' pN ms/um ']);
[k1, dk1, sig]=round_significance(gamma_psd_II*1e9, sigma_gamma_psd_II*1e9);
disp(['gamma_psd:' k1 '\pm' dk1])

disp(['D_psd: ' num2str(D_psd_II*1e12) '+-' num2str(sigma_D_psd_II*1e12) ' um^2/s'])
[k1, dk1, sig]=round_significance(D_psd_II*1e12, sigma_D_psd_II*1e12);
disp(['D_psd: ' k1 '\pm' dk1])
disp('................')
%%
titleIII='Experiment III, P=9.2mW';
[k_psd_III, sigma_k_psd_III, gamma_psd_III, sigma_gamma_psd_III,  D_psd_III, sigma_D_psd_III]=plotsub_psd('Data_positions_Fig9_1P6_S.mat',[bx1+2*bx2+2*xwi 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix, titleIII, T, subs,1);


disp('................')
disp(titleIII)
disp('Power spectrum density analysis by linear fitting')
disp(['k_psd: ' num2str(k_psd_III*1e6) '+-' num2str(sigma_k_psd_III*1e6) 'pN/um'])
[k1, dk1, sig]=round_significance(k_psd_III*1e6, sigma_k_psd_III*1e6);
disp(['k_psd: ' k1 '\pm' dk1 ])

disp(['gamma_psd:' num2str(gamma_psd_III*1e9) '+-'  num2str(sigma_gamma_psd_III*1e9) ' pN ms/um ']);
[k1, dk1, sig]=round_significance(gamma_psd_III*1e9, sigma_gamma_psd_III*1e9);
disp(['gamma_psd:' k1 '\pm' dk1])

disp(['D_psd: ' num2str(D_psd_III*1e12) '+-' num2str(sigma_D_psd_III*1e12) ' um^2/s'])
[k1, dk1, sig]=round_significance(D_psd_III*1e12, sigma_D_psd_III*1e12);
disp(['D_psd: ' k1 '\pm' dk1])
disp('................')
