close all, clear all;

%%  ==============Parameter declaration============


kB=1.38e-23; % Boltzmann constant [m^2kg/s^2K]
T=300;  % Temperature [K]
r=1.03E-6;      % Particle radius [m]
v=0.00002414*10^(247.8/(-140+T));  % Water viscosity [Pa*s]
gamma=pi*6*r*v; %[m*Pa*s]



%%  =========Loading selected file============
load('Data_positions_Fig9_1P4_S.mat')

addpath pot


%Boltzmann constant
kb=1.38064852e-23;

%number of bins of the histogram, if not set default is 50
P=100; 

%use a subsampled data set
subs=1;

%linear fit
[k_pot_lf, sigma2_k_pot_lf, x_alpha_lf, mrho_lf, sigma2_rho_lf, mU_lf, sigma2_Ulf, rho0_lf, x_eq_lf, U_0_lf]=pot_lfit(x(1:subs:size(x,1),:),T,P);

%non-linear fit
[k_pot_nl, sigma2_k_potmnl, x_alpha_nl, mrho_nl, sigma2_rho_nl, mU_nl, sigma2_U_nl, rho0_nl, x_eq_nl]=pot_nlfit(x(1:subs:size(x,1),:),T,P);

%% plot figures

col1=[0,0.6,0.8];

col2=[0.9,0.4,0.2];
xwi = 400;    % width of the plot square
bx1 = 80;     % extra space at the left
bx2 = 20;     % extra space at the right

Xpix = 2*xwi+2*bx1+2*bx2;  % total

ywi = 300;    % length riquadro con funzione
by1 = 60;     % extra space below
by2 = 30;     % extra space up

Ypix = by1+ywi+by2;  % larghezza figura in pixel


figure('Position',[10 20 Xpix Ypix]); % crea la figura
axes( 'Position',[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  % fa in modo di centrare il riquadro degli assi nella posizione voluta


bar(x_alpha_lf, mrho_lf)
hold on 
rhomodel_lf=rho0_lf*exp(-k_pot_lf/(2*kb*T)*(x_alpha_lf-x_eq_lf).^2);
rhomodel_nl=rho0_nl*exp(-k_pot_nl/(2*kb*T)*(x_alpha_nl-x_eq_nl).^2);

plot(x_alpha_lf, rhomodel_lf,'LineWidth',3,'Color',col2)
plot(x_alpha_nl, rhomodel_nl,'LineWidth',3,'Color','green')
box on
xticks((-0.5:0.1:0.5)*1e-7);
xlim([x_alpha_lf(1) x_alpha_lf(end)]);

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',15);
xlabel('$x(\mu m)$','Interpreter','Latex', 'FontSize',20)
ylabel('$\rho (\rm {counts})$','Interpreter','Latex', 'FontSize',20)
  


axes('Position',[(bx2+2*bx1+xwi) 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
hold on;

scatter(x_alpha_lf, mU_lf/(kB*T),60,'o')
U_model_lf=-(log(rhomodel_lf));
U_model_nl=-(log(rhomodel_nl));
scatter(x_alpha_lf,U_model_lf,60,'o','markerfacecolor', [0.00,0.45,0.74],'markeredgecolor',[0.00,0.45,0.74])



hold on

plot(x_alpha_lf,U_model_lf, 'LineWidth',3,'Color',col2)
plot(x_alpha_nl,U_model_nl, 'LineWidth',3,'Color','green')
box on

xticks((-0.5:0.1:0.5)*1e-7);
xlim([x_alpha_lf(1) x_alpha_lf(end)]);
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5, 'FontSize',15);
xlabel('$x(\mu m)$','Interpreter','Latex', 'FontSize',20)
ylabel('$U(k_BT)$','Interpreter','Latex','FontSize',20)

  