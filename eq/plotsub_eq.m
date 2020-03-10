function [k_pot_lf, k_pot_nl]=plotsub_eq(filename, positioninthefig1, positionintefig2, title1, T, P)
load(filename);
disp(filename);
addpath pot
kb = 1.38064852e-23;

%translate everithing to zero
x = x - repmat(mean(x),size(x,1),1);

%default number of bins
P=50;

%user defined number of bins


[x_alpha, mrho, sigma2_rho, frequencynor, logf, mlogf, sigma_logh,mU,  sigma2_U, mloghist, Eloghist, U_0_exp]=prob_dist_energy(x,P, T)

%delete zeros to avoid Inf in weights
sigma2_rho(sigma2_rho==0) = 1;







% Guess for the initial conditions for the non-linear fitting
a0=(max(mrho)*sqrt(pi))^2;
b0=0;



%blue color
col1=[73/255,4/255,10/255];
	
%yellow
col2=[241/255,185/255,14/255];
%gray color for experimental data
	

%colbar=[4/255,45/255,73/255];
	colbar=[7/255, 79/255, 129/255];
% xwi = 400;    % width of the plot square
% bx1 = 80;     % extra space at the left
% bx2 = 20;     % extra space at the right
% 
% Xpix = 3.5*xwi+3*bx1+3*bx2;  % total
% 
% ywi = 300;    % length riquadro con funzione
% by1 = 60;     % extra space below
% by2 = 30;     % extra space up
% Ypix = 2*by1+2*ywi+3*by2;  % larghezza figura in pixel

P=50;
subsample=1;
%number of bins of the histogram, if not set default is 50
%linear fit
[k_eq,sigma2_k_eq]=eq1d(x(1:subsample:1000000,:),T,0e-9);
axes( 'Position',positioninthefig1);  % fa in modo di centrare il riquadro degli assi nella posizione voluta

rhomodel=a0*exp(-k_eq/(2*kb*T)*(x_alpha).^2);


plot(x_alpha*1e6, rhomodel*1e-6, 'LineWidth',3,'Color',col2, 'DisplayName',  'Equipotential');
hold on 
errorbar(x_alpha*1e6,  mrho*1e-6, 1e-6*abs(sigma2_rho),'o','MarkerSize',7 ,'LineWidth', 1.5,'Color',colbar, 'DisplayName', 'Experimental probability distribution');
box on

%xticks((-0.5:0.1:0.5)*1e-7);
xlim([-0.06 0.06]);
ylim([0, 47.5])
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',15);
xlabel('$x(\mu m)$','Interpreter','Latex', 'FontSize',20)
ylabel('$\rho (\rm{\mu m^{-1}})$','Interpreter','Latex', 'FontSize',20)
hold off
title(title1)
legend

%second figure, Energy potential distribution, exp I
axes('Position',positionintefig2);  
hold on;

%scatter(x_alpha_lf*1e6, -log(mrho_lf)-U_0_exp,80,'o', 'markerfacecolor','colbar', 'markeredgecolor', colbar , 'DisplayName', 'Experimental values of potential energy')
%errorbar(x_alpha_lf*1e6,  -log(mrho_lf)-U_0_exp_nl, -log(1e-6*abs(sigma2_rho_lf)), 'Color', colbar,  'DisplayName', 'Experimental values of potential energy');
%sigma2_U_nl


U_model_lf=-(log(rhomodel));





plot(x_alpha*1e6,U_model_lf-U_0_exp, 'LineWidth',3,'Color',col2,'DisplayName',  'Linear fitting')
hold on



errorbar(x_alpha*1e6,  -log(mrho)-U_0_exp, sigma2_U/(kb*T), 'o','MarkerSize',7,'LineWidth', 1.5, 'Color', colbar,  'DisplayName', 'Experimental values of potential energy');

xlim([-0.06 0.06]);
%xlim([x_alpha_lf(1)*1e6 x_alpha_lf(end)*1e6]);
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5, 'FontSize',15);
xlabel('$x(\mu m)$','Interpreter','Latex', 'FontSize',20)
ylabel('$U(k_BT)$','Interpreter','Latex','FontSize',20)
ylim([0, 10])

hold off

legend
end