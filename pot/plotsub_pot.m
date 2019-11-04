function [k_pot_lf,sigma2_k_pot_lf, k_pot_nl,sigma2_k_pot_nl , k_eq, sigma2_k_eq]=plotsub_pot(filename, positioninthefig1, positionintefig2, nbins_pot, subs)
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

%number of bins of the histogram, if not set default is 50
%linear fit
[k_pot_lf, sigma2_k_pot_lf, x_alpha_lf, mrho_lf, sigma2_rho_lf, mU_lf, sigma2_U_lf, rho0_lf, x_eq_lf, U_0_exp_lf]=pot_lfit(x(1:subs:end,:),T,nbins_pot);
%non-linear fit`
[k_pot_nl, sigma2_k_pot_nl, x_alpha_nl, mrho_nl, sigma2_rho_nl, mU_nl, sigma2_U_nl, rho0_nl, x_eq_nl,  U_0_exp_nl]=pot_nlfit(x(1:subs:end,:),T,nbins_pot);

 
%number of bins of the histogram, if not set default is 50
%linear fit
[k_eq,sigma2_k_eq]=eq1d(x(1:subs:end,:),T);
 
axes( 'Position',positioninthefig1);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
% disp('Check that the probability distribiution is normalized')
% disp('Integral of the experimental distribuitions')
% disp('Linear fitting')
% disp(sum(mrho_lf*1e-6)*((x_alpha_lf(2)-x_alpha_lf(1))*1e6));
% disp('Non-linear fitting')
% disp(sum(mrho_lf*1e-6)*((x_alpha_lf(2)-x_alpha_lf(1))*1e6));
 
a0=(max(mrho_lf))^2;
rhomodel_eq=rho0_lf*exp(-k_eq/(2*kb*T)*(x_alpha_lf).^2);
rhomodel_lf=rho0_lf*exp(-k_pot_lf/(2*kb*T)*(x_alpha_lf-x_eq_lf).^2);
rhomodel_nl=rho0_nl*exp(-k_pot_nl/(2*kb*T)*(x_alpha_nl-x_eq_nl).^2);
% disp('Integral of the fitted probability distributions')
% disp('Linear fitting')
% disp(sum(rhomodel_lf*1e-6)*((x_alpha_lf(2)-x_alpha_lf(1))*1e6));
% disp('Non-linear fitting')
% disp(sum(rhomodel_nl*1e-6)*((x_alpha_lf(2)-x_alpha_lf(1))*1e6));
cyan=[0/255,128/255,128/255];
scatter((x_alpha_lf)*1e6+0.0005, rhomodel_eq*1e-6,'SizeData',200,'MarkerEdgeColor','white', 'MarkerFaceColor',color2rgb('white_cyan'),'MarkerFaceAlpha',0.9,'DisplayName',  'Equipartition');
hold on
 
plot(x_alpha_lf*1e6, rhomodel_lf*1e-6, 'LineWidth',4.5,'Color',col2, 'DisplayName',  'Linear fit');
 
cyan=[0/255,128/255,128/255];
plot(x_alpha_nl*1e6, rhomodel_nl*1e-6,'--','LineWidth',3,'Color',color2rgb('deep_purple'), 'DisplayName', 'Non-linear fit')
 
errorbar(x_alpha_lf*1e6,  mrho_lf*1e-6, 1e-6*abs(sigma2_rho_lf),'.','MarkerSize',7 ,'LineWidth', 1.5,'Color',colbar,'DisplayName', 'Experimental');
 
 
 
box on
%xticks((-0.5:0.1:0.5)*1e-7);
xlim([-0.070 0.070]);
ylim([-1 50])
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5,'FontSize',25);
xlabel(' x( Âµ m)','Interpreter','tex', 'FontSize',30,'FontName','Times New Roman');
 
%xlabel('$x (Âµ m)$','Interpreter','tex', 'FontSize',30)
%xlabel(['c_{O2} [',char(181),'mol/l]'],'Interpreter','tex')
ylabel(['\rho (Âµ m^{-1})'],'Interpreter','tex', 'FontSize',30, 'FontName','Times New Roman')
hold off
legend('FontSize',25)
 
%second figure, Energy potential distribution, exp I
axes('Position',positionintefig2);  
hold on;
 

 
 
U_model_lf=-(log(rhomodel_lf));
U_model_nl=-(log(rhomodel_nl));
 
 
U_model_eq=-(log(rhomodel_eq));
 
scatter((x_alpha_lf)*1e6+0.0005,U_model_eq-U_0_exp_lf, 'SizeData',200,'MarkerEdgeColor','white','MarkerFaceColor',color2rgb('white_cyan'),'MarkerFaceAlpha',0.9, 'DisplayName',  'Equipartition')
hold on
plot(x_alpha_lf*1e6,U_model_lf-U_0_exp_lf, 'LineWidth',4.5,'Color',col2,'DisplayName',  'Linear fit')
 
plot(x_alpha_nl*1e6,U_model_nl-U_0_exp_nl, '--', 'LineWidth',3,'Color',color2rgb('deep_purple'), 'DisplayName', 'Non-linear fit')
 
%plot(x_alpha_lf*1e6,U_model_eq-U_0_exp_lf,'-.', 'LineWidth',2,'Color','black','DisplayName',  'Equipotential')
 
box on
 
%xticks((-0.5:0.1:0.5)*1e-7);
 
errorbar(x_alpha_lf*1e6,  -log(mrho_lf)-U_0_exp_nl, sigma2_U_lf/(kb*T), '.','MarkerSize',7,'LineWidth', 1.5, 'Color', colbar,  'DisplayName', 'Experimental values of potential energy');
 
xlim([-0.070 0.070]);
%xlim([x_alpha_lf(1)*1e6 x_alpha_lf(end)*1e6]);
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
 
xlabel(' x( Âµ m)','FontSize',30,'FontName','Times New Roman');
% 
% 
 
% 
ylabel('$U(k_{\rm B} T)$','Interpreter','Latex','FontSize',30)
ylim([-0.5, 10])
 
hold off
 set(gca,'DefaultTextFontname', 'Times')
   set(gca,'DefaultAxesFontName','Times')
   legend('FontSize',25)
end
