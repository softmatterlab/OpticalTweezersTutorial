% Potential analysis for Nexp

% Initialization of the workspace
clear;

close all;

addpath pot

% Load data file

load('Data_positions_Fig9_1P4_S.mat');

%Vx=Vx*Sdx;

kb=1.38064852e-23;

subs=1; %use a subsampled data set

%linear fit
[xbins,mU,EU,k_pot,Ek_pot,mhist,Ehist,h0,x_eq]=pot_lfit_v1(x(1:subs:size(x,1),1),T,50);
disp(['tau_0: ' num2str(6*pi*eta*a/k_pot)])
disp(['dt: ' num2str(dt*subs)])
%non-linear fit
[xbins,mU,EU,k_pot,Ek_pot,mhist,Ehist,h0,x_eq]=pot_nlfit_v1(x(1:subs:size(x,1),1),T,50);


disp(['tau_0: ' num2str(6*pi*eta*a/k_pot)])
disp(['dt: ' num2str(dt*subs)])


