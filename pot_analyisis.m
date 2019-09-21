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
[k_potlf, Ek_potlf, xbinslf, mrholf, Erholf, mUlf, EUlf, rho0lf, x_eqlf]=pot_lfit(x(1:subs:size(x,1),:),T,100);

%non-linear fit
[k_potnl, Ek_potml, xbinsnl, mrhonl, Erhonl, mUnl, EUnl, rho0nl, x_eqnl]=pot_nlfit(x(1:subs:size(x,1),:),T,100);

