close all, clear all;

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

Ypix = 2*by1+2*ywi+3*by2;  % larghezza figura in pixel
%number of bins of the histogram, if not set default is 50
P=50; 

%use a subsampled data set
subs=1;

%%  =========Loading selected file============
%load('Data_positions_Fig9_1P6_S.mat')

addpath pot


%Boltzmann constant
kb=1.38064852e-23;




%% plot figures



%creates the figure to do the subplots
figure('Position',[10 20 Xpix Ypix]);
%first figure, probability distribution, exp I

%%
titleI='Experiment I, P=2.3mW';
[k_pot_lf_I, k_pot_nl_I]=plotsub('Data_positions_Fig9_1P2_S.mat',[bx1 0 xwi 0]/Xpix + [0 2*by1+ywi+by2 0 ywi]/Ypix, [bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix, titleI, T, P, subs)

titleII='Experiment II, P=6.0mW';
[k_pot_lf_II, k_pot_nl_II]=plotsub('Data_positions_Fig9_1P4_S.mat',[2*bx1+xwi+bx2 0 xwi 0]/Xpix + [0 2*by1+ywi+by2 0 ywi]/Ypix, [2*bx1+xwi+bx2 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix, titleII, T, P, subs)


titleIII='Experiment III, P=9.2mW';
[k_pot_lf_III, k_pot_nl_III]=plotsub('Data_positions_Fig9_1P6_S.mat',[3*bx1+2*xwi+2*bx2 0 xwi 0]/Xpix + [0 2*by1+ywi+by2 0 ywi]/Ypix, [3*bx1+2*xwi+2*bx2 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix, titleIII, T, P, subs)

%2:3mW, 6:0mW, and 9:2mW