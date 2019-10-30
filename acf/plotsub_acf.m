function [k_acf_lf,Ek_acf_lf,D_acf_lf,ED_acf_lf,gamma_acf_lf, sigma2_gamma_acf_lf , k_acf_nl, Ek_acf_nl, D_acf_nl, ED_acf_nl,gamma_acf_nl, sigma2_gamma_acf_nl, tau0_exp_lf, tau0_exp_nl]=plotsub_pot(filename, positioninthefig1, title1, T, subsample)
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
[k_acf_lf, Ek_acf_lf, D_acf_lf, ED_acf_lf,gamma_acf_lf, sigma2_gamma_acf_lf,tau_acf_lf, mc, Ec,indc, tau0_exp_lf, c0_exp_lf]=acf_lfit(x(1:subsample:size(x,1),:),T,dt*subsample);
%[k_acf_lf, Ek_acf, D_acf, ED_acf, tau_lf, mc, Ec,indc, tau0_exp_lf, c0_exp_lf]=acf_lfit(x(1:subsample:size(x,1),:),T,dt*subsample);
%non linear fit
[k_acf_nl, Ek_acf_nl, D_acf_nl, ED_acf_nl,gamma_acf_nl, sigma2_gamma_acf_nl, tau_nl, mc, Ec, indc, tau0_exp_nl, c0_exp_nl]=acf_nlfit(x(1:subsample:size(x,1),:),T,dt*subsample);

axes( 'Position',positioninthefig1);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
%bar(x_alpha_lf*1e6, mrho_lf*1e-6, 'EdgeColor', 'white', 'FaceColor', colbar, 'facealpha', 0.8, 'HandleVisibility','off')




plot(tau_acf_lf(1:20:15*indc),c0_exp_lf*exp(-tau_acf_lf(1:20:15*indc)/tau0_exp_lf)*1e12, 'LineWidth',3,'Color',col2, 'DisplayName',  'Linear fitting')
hold on
plot(tau_nl(1:20:15*indc),c0_exp_nl*exp(-tau_nl(1:20:15*indc)/tau0_exp_nl)*1e12, '--', 'LineWidth',3,'Color',color2rgb('deep_purple'), 'DisplayName',  'Non -linear fitting')
errorbar(tau_acf_lf(1:20:15*indc),mc(1:20:15*indc)*1e12,Ec(1:20:15*indc)*1e12,'.','MarkerSize',1 ,'LineWidth', 1.5,'Color',colbar, 'DisplayName', 'Experimental autocorrelation function');
% 

%plot(tau(1:20:6*indc),c0_exp*exp(tau(1:20:6*indc)/tau0_exp)*1e12,'b')
% % 
% disp('Integral of the fitted probability distributions')
% disp('Linear fitting')
% disp(sum(rhomodel_lf*1e-6)*((x_alpha_lf(2)-x_alpha_lf(1))*1e6));
% disp('Non-linear fitting')
% disp(sum(rhomodel_nl*1e-6)*((x_alpha_lf(2)-x_alpha_lf(1))*1e6));
% plot(x_alpha_lf*1e6, rhomodel_lf*1e-6, 'LineWidth',3,'Color',col2, 'DisplayName',  'Linear fitting');
% hold on 
% 
% plot(x_alpha_nl*1e6, rhomodel_nl*1e-6,'--','LineWidth',3,'Color',col1, 'DisplayName', 'Non-linear fitting')
% errorbar(x_alpha_lf*1e6,  mrho_lf*1e-6, 1e-6*abs(sigma2_rho_lf),'o','MarkerSize',7 ,'LineWidth', 1.5,'Color',colbar, 'DisplayName', 'Experimental probability distribution');
% box on
% %xticks((-0.5:0.1:0.5)*1e-7);
ntaus=6;
plot([tau0_exp_lf*ntaus,tau0_exp_lf*ntaus],[-0.5,1],'--k', 'HandleVisibility','off')
text(tau0_exp_lf*ntaus*0.95,0.3e-4,[num2str(ntaus),'$\tau_0$'],'Interpreter','latex','FontSize',30)
xlim([0 0.009]);
ylim([-0.1, 3.1]*1e-4)
% set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',15);
% xlabel('$x(\mu m)$','Interpreter','Latex', 'FontSize',20)
% ylabel('$\rho (\rm{\mu m^{-1}})$','Interpreter','Latex', 'FontSize',20)
% hold off
 
% legend
% 
% %second figure, Energy potential distribution, exp I
% axes('Position',positionintefig2);  
% hold on;
% 
% %scatter(x_alpha_lf*1e6, -log(mrho_lf)-U_0_exp,80,'o', 'markerfacecolor','colbar', 'markeredgecolor', colbar , 'DisplayName', 'Experimental values of potential energy')
% %errorbar(x_alpha_lf*1e6,  -log(mrho_lf)-U_0_exp_nl, -log(1e-6*abs(sigma2_rho_lf)), 'Color', colbar,  'DisplayName', 'Experimental values of potential energy');

% 
% 
% 
% plot(x_alpha_lf*1e6,U_model_lf-U_0_exp_lf, 'LineWidth',3,'Color',col2,'DisplayName',  'Linear fitting')
% hold on
% plot(x_alpha_nl*1e6,U_model_nl-U_0_exp_nl, '--', 'LineWidth',3,'Color',col1, 'DisplayName', 'Non-linear fitting')
% box on
% 
% %xticks((-0.5:0.1:0.5)*1e-7);
% 
% errorbar(x_alpha_lf*1e6,  -log(mrho_lf)-U_0_exp_nl, sigma2_U_lf/(kb*T), 'o','MarkerSize',7,'LineWidth', 1.5, 'Color', colbar,  'DisplayName', 'Experimental values of potential energy');
% 
% xlim([-0.06 0.06]);
% %xlim([x_alpha_lf(1)*1e6 x_alpha_lf(end)*1e6]);
 set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5, 'FontSize',25);
 xlabel('$\tau(\rm s)$','Interpreter','Latex', 'FontSize',30)
 ylabel('$C(  \mu  m ^2)$','Interpreter','Latex','FontSize',30, 'FontName', 'TimesNewRoman')
% ylim([0, 10])
% 
% hold off
% 
 legend
end