function [k_pot_lf, k_pot_nl]=plotsub(filename, positioninthefig1, positionintefig2, title1, T, P, subssample)
load(filename)
kb=1.38e-23;
col1=[232/255,193/255,95/255];
%blue color
col2=[33/255,132/255,218/255];
%gray color for experimental data
colbar=[119/255,136/255,153/255];
xwi = 400;    % width of the plot square
bx1 = 80;     % extra space at the left
bx2 = 20;     % extra space at the right

Xpix = 3.5*xwi+3*bx1+3*bx2;  % total

ywi = 300;    % length riquadro con funzione
by1 = 60;     % extra space below
by2 = 30;     % extra space up


P=50;
Ypix = 2*by1+2*ywi+3*by2;  % larghezza figura in pixel
%number of bins of the histogram, if not set default is 50
%linear fit
[k_pot_lf, sigma2_k_pot_lf, x_alpha_lf, mrho_lf, sigma2_rho_lf, mU_lf, sigma2_Ulf, rho0_lf, x_eq_lf, U_0_lf]=pot_lfit(x(1:subssample:size(x,1),:),T,P);

%non-linear fit
[k_pot_nl, sigma2_k_potmnl, x_alpha_nl, mrho_nl, sigma2_rho_nl, mU_nl, sigma2_U_nl, rho0_nl, x_eq_nl]=pot_nlfit(x(1:subssample:size(x,1),:),T,P);

axes( 'Position',positioninthefig1);  % fa in modo di centrare il riquadro degli assi nella posizione voluta


bar(x_alpha_lf*1e6, mrho_lf, 'EdgeColor', 'white', 'FaceColor', colbar, 'facealpha', 0.5, 'DisplayName', 'Experimental probability distribution')
hold on 
rhomodel_lf=rho0_lf*exp(-k_pot_lf/(2*kb*T)*(x_alpha_lf-x_eq_lf).^2);
rhomodel_nl=rho0_nl*exp(-k_pot_nl/(2*kb*T)*(x_alpha_nl-x_eq_nl).^2);

plot(x_alpha_lf*1e6, rhomodel_lf,'LineWidth',3,'Color',col2, 'DisplayName',  'Linear fitting');
plot(x_alpha_nl*1e6, rhomodel_nl,'LineWidth',3,'Color',col1, 'DisplayName', 'Non-linear fitting')
box on
%xticks((-0.5:0.1:0.5)*1e-7);
xlim([x_alpha_lf(1)*1e6 x_alpha_lf(end)*1e6]);

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',15);
xlabel('$x(\mu m)$','Interpreter','Latex', 'FontSize',20)
ylabel('$\rho (\rm {counts})$','Interpreter','Latex', 'FontSize',20)
hold off
title(title1)
legend

%second figure, Energy potential distribution, exp I
axes('Position',positionintefig2);  
hold on;

scatter(x_alpha_lf*1e6, -log(mrho_lf),60,'o', 'markerfacecolor', colbar, 'markeredgecolor', colbar , 'DisplayName', 'Experimental values of potential energy')
U_model_lf=-(log(rhomodel_lf));
U_model_nl=-(log(rhomodel_nl));


hold on

plot(x_alpha_lf*1e6,U_model_lf, 'LineWidth',3,'Color',col2,'DisplayName',  'Linear fitting')
plot(x_alpha_nl*1e6,U_model_nl, 'LineWidth',3,'Color',col1, 'DisplayName', 'Non-linear fitting')
box on

%xticks((-0.5:0.1:0.5)*1e-7);
xlim([x_alpha_lf(1)*1e6 x_alpha_lf(end)*1e6]);
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5, 'FontSize',15);
xlabel('$x(\mu m)$','Interpreter','Latex', 'FontSize',20)
ylabel('$U(k_BT)$','Interpreter','Latex','FontSize',20)


hold off

legend
end