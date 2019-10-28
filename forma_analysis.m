clear all 
close all
%load data files
load('Data_positions_Fig9_1P2_S.mat')

xx = reshape(x,[size(x,1)*size(x,2),1]);
N=length(xx);
%some useful calculations from the positions


%physical constants
kB=1.38e-23; % Boltzmann constant [m^2kg/s^2K]

%definition of parameter

eta=0.00002414*10^(247.8/(-140+T));  % Water viscosity [Pa*s]
k_th=14*1e-6;
gamma=6*pi*eta*a;

addpath forma


[k_forma1d, D_forma1d] = forma1d(xx, dt, gamma);

disp(['FORMA 1D results:'])
disp(['k* = ' num2str( k_forma1d ) ' N/m'])
disp(['k*/k = ' num2str( k_forma1d/k )])
disp(['D* = ' num2str( D_forma1d ) ' m^2/s'])
disp(['D*/D = ' num2str( D_forma1d/D )])