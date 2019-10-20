function [k_pot_lf, k_pot_nl]=plotsub_msd(filename, positioninthefig1, positionintefig2, title1, T, P, subssample)
load(filename);
disp(filename);
kb=1.38e-23;

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

%number of bins of the histogram, if not set default is 50
%linear fit
%[k_pot_lf, ~, x_alpha_lf, mrho_lf, sigma2_rho_lf,~,~, rho0_lf, x_eq_lf,U_0_exp]=pot_lfit(x(1:subssample:size(x,1),:),T,P);
[k_pot_lf, sigma2_k_pot_lf, x_alpha_lf, mrho_lf, sigma2_rho_lf, mU_lf, sigma2_U_lf, rho0_lf, x_eq_lf, U_0_exp_lf]=pot_lfit(x,T,P);
%non-linear fit`
%[k_pot_nl, ~, x_alpha_nl, ~, ~, ~, ~, rho0_nl, x_eq_nl,~]=pot_nlfit(x(1:subssample:size(x,1),:),T,P);
[k_pot_nl, sigma2_k_pot_nl, x_alpha_nl, mrho_nl, sigma2_rho_nl, mU_nl, sigma2_U_nl, rho0_nl, x_eq_nl,  U_0_exp_nl]=pot_nlfit(x,T,P)
axes( 'Position',positioninthefig1);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
disp('Check that the probability distribiution is normalized')
disp('Integral of the experimental distribuitions')
disp('Linear fitting')
disp(sum(mrho_lf*1e-6)*((x_alpha_lf(2)-x_alpha_lf(1))*1e6));
disp('Non-linear fitting')
disp(sum(mrho_lf*1e-6)*((x_alpha_lf(2)-x_alpha_lf(1))*1e6));
%bar(x_alpha_lf*1e6, mrho_lf*1e-6, 'EdgeColor', 'white', 'FaceColor', colbar, 'facealpha', 0.8, 'HandleVisibility','off')


rhomodel_lf=rho0_lf*exp(-k_pot_lf/(2*kb*T)*(x_alpha_lf-x_eq_lf).^2);
rhomodel_nl=rho0_nl*exp(-k_pot_nl/(2*kb*T)*(x_alpha_nl-x_eq_nl).^2);
disp('Integral of the fitted probability distributions')
disp('Linear fitting')
disp(sum(rhomodel_lf*1e-6)*((x_alpha_lf(2)-x_alpha_lf(1))*1e6));
disp('Non-linear fitting')
disp(sum(rhomodel_nl*1e-6)*((x_alpha_lf(2)-x_alpha_lf(1))*1e6));
plot(x_alpha_lf*1e6, rhomodel_lf*1e-6, 'LineWidth',3,'Color',col2, 'DisplayName',  'Linear fitting');
hold on 

plot(x_alpha_nl*1e6, rhomodel_nl*1e-6,'--','LineWidth',3,'Color',col1, 'DisplayName', 'Non-linear fitting')
errorbar(x_alpha_lf*1e6,  mrho_lf*1e-6, 1e-6*abs(sigma2_rho_lf),'o','MarkerSize',7 ,'LineWidth', 1.5,'Color',colbar, 'DisplayName', 'Experimental probability distribution');
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


U_model_lf=-(log(rhomodel_lf));
U_model_nl=-(log(rhomodel_nl));




plot(x_alpha_lf*1e6,U_model_lf-U_0_exp_lf, 'LineWidth',3,'Color',col2,'DisplayName',  'Linear fitting')
hold on
plot(x_alpha_nl*1e6,U_model_nl-U_0_exp_nl, '--', 'LineWidth',3,'Color',col1, 'DisplayName', 'Non-linear fitting')
box on

%xticks((-0.5:0.1:0.5)*1e-7);

errorbar(x_alpha_lf*1e6,  -log(mrho_lf)-U_0_exp_nl, sigma2_U_lf/(kb*T), 'o','MarkerSize',7,'LineWidth', 1.5, 'Color', colbar,  'DisplayName', 'Experimental values of potential energy');

xlim([-0.06 0.06]);
%xlim([x_alpha_lf(1)*1e6 x_alpha_lf(end)*1e6]);
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5, 'FontSize',15);
xlabel('$x(\mu m)$','Interpreter','Latex', 'FontSize',20)
ylabel('$U(k_BT)$','Interpreter','Latex','FontSize',20)
ylim([0, 10])

hold off

legend
end