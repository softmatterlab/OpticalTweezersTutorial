close all;
clear all;

%%  ==============Parameter declaration============


kB=1.38e-23; % Boltzmann constant [m^2kg/s^2K]
T=300;  % Temperature [K]
r=1.03E-6;      % Particle radius [m]
v=0.00002414*10^(247.8/(-140+T));  % Water viscosity [Pa*s]
gamma=pi*6*r*v; %[m*Pa*s]

xwi = 400;    % width of the plot square
bx1 = 80;     % extra space at the left
bx2 = 20;     % extra space at the right

Xpix = 3*xwi+3*bx1+3*bx2;  % total

ywi = 300;    % length riquadro con funzione
by1 = 60;     % extra space below
by2 = 30;     % extra space up

Ypix = 1*by1+1*ywi+1*by2;  % larghezza figura in pixel


%%  =========Loading selected file============
%load('Data_positions_Fig9_1P6_S.mat')

addpath msd
addpath wlsice


%Boltzmann constant
kb=1.38064852e-23;




%% plot figures



%creates the figure to do the subplots
figure('Position',[10 20 Xpix Ypix]);
%first figure, probability distribution, exp I

%%

%subs=10; %use a subsampled data set
%exp I SUBS=20, maxlag=50 
%exp II subs=17, maxlag=25
%exp III subs=10, maxlag=25
titleI='Experiment I, P=2.3mW';
[k_msd_I]=plotsub_msd('Data_positions_Fig9_1P2_S.mat',[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix,  titleI, T, 20, 50);
%%
titleII='Experiment II, P=6.0mW';
[k_msd_II]=plotsub_msd('Data_positions_Fig9_1P4_S.mat',[2*bx1+xwi+bx2 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix, titleII, T, 17, 70);


titleIII='Experiment III, P=9.2mW';
[k_msd_III]=plotsub_msd('Data_positions_Fig9_1P6_S.mat',[3*bx1+2*xwi+2*bx2 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix, titleIII, T, 10, 100);

%2:3mW, 6:0mW, and 9:2mW