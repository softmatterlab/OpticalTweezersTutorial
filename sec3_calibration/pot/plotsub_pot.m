function [k_pot_lf,sigma2_k_pot_lf, k_pot_nl,sigma2_k_pot_nl , k_eq, sigma2_k_eq]=plotsub_pot(filename, positioninthefig1, positionintefig2, title1, T, nbins_pot, subssample,aa)
%function [k_pot_lf,sigma2_k_pot_lf, k_pot_nl,sigma2_k_pot_nl , k_eq, sigma2_k_eq]=plotsub_pot(filename, positioninthefig1, positionintefig2, title1, T, nbins_pot, subssample,aa)
%computes the potential and equipartition method analysis and plots the fit
%it calls the functions pot_lfit, pot_nlfit and eq1d
load(filename);
disp(filename);
kb=1.38064852e-23;

col1=[73/255,4/255,10/255];

col2=[241/255,185/255,14/255];

col2=[1,0,0];

col3=[0.00,0.45,0.74];

colbar=[7/255, 79/255, 129/255];


%number of bins of the histogram, if not set default is 50
%linear fit
[k_pot_lf, sigma2_k_pot_lf, x_alpha_lf, mrho_lf, sigma2_rho_lf, mU_lf, sigma2_U_lf, rho0_lf, x_eq_lf, U_0_exp_lf]=pot_lfit(x(1:subssample:end,:),T,nbins_pot);
%non-linear fit`
[k_pot_nl, sigma2_k_pot_nl, x_alpha_nl, mrho_nl, sigma2_rho_nl, mU_nl, sigma2_U_nl, rho0_nl, x_eq_nl,  U_0_exp_nl]=pot_nlfit(x(1:subssample:end,:),T,nbins_pot);

 
%number of bins of the histogram, if not set default is 50
%linear fit
[k_eq,sigma2_k_eq]=eq1d(x(1:subssample:end,:),T, 0);
 
axes( 'Position',positioninthefig1);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
hold on
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

e=errorbar(x_alpha_lf*1e9,  mrho_lf*1e-9, 1e-9*abs(sigma2_rho_lf),'.','MarkerSize',30 ,'LineWidth', 1.5,'Color',colbar,'DisplayName', 'Experimental');
e.Color = col3;
hold on
scatter((x_alpha_lf)*1e9+0.0005, rhomodel_eq*1e-9,155,'o','markeredgecolor','k','DisplayName', 'Equipartition');


% scatter(tau_acf_lf(1:20:end_plot*indc),mc(1:20:end_plot*indc)*1e12,60,'o','markerfacecolor', col3,'markeredgecolor',col3);


plot(x_alpha_lf*1e9, rhomodel_lf*1e-9, 'LineWidth',3,'Color','k', 'DisplayName',  'Linear fit');



cyan=[0/255,128/255,128/255];
plot(x_alpha_nl(2:end-1)*1e9, rhomodel_nl(2:end-1)*1e-9,'--','LineWidth',3,'Color','r', 'DisplayName',  'Non -linear fitting')

box on

xlim([-70 70]);
ylim([-1e-3 60e-3])

if aa==5
    set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5,'FontSize',25);
else
    set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5,'FontSize',25,'Yticklabel',[]);
end

set(gca,'Xticklabel',[]);

if aa==5
    ylabel('$\rho ( \rm nm^{-1}$)','Interpreter','Latex', 'FontSize',30);
end

if aa==5
    LL= legend({'Experimental', 'Equipartition','Linear fitting', 'Non-linear fitting'},'Box','off','Position',[0.2 0.73 0.1 0.2]);
    LL.FontSize = 18;
    
end

hold off
axes('Position',positionintefig2);
hold on;
 

 
 
U_model_lf=-(log(rhomodel_lf));
U_model_nl=-(log(rhomodel_nl));
 
 
U_model_eq=-(log(rhomodel_eq));

e=errorbar(x_alpha_lf*1e9,  -log(mrho_lf)-U_0_exp_nl, sigma2_U_lf/(kb*T), '.','MarkerSize',30,'LineWidth', 1.5, 'Color', colbar,  'DisplayName', 'Experimental values of potential energy');
e.Color = col3;
xlim([-70 70]);

scatter((x_alpha_lf)*1e9+0.0005,U_model_eq-U_0_exp_lf, 155,'o','markeredgecolor','k', 'DisplayName',  'Equipartition')



hold on
plot(x_alpha_lf*1e9,U_model_lf-U_0_exp_lf, 'LineWidth',3,'Color','k', 'DisplayName',  'Linear fit');


plot(x_alpha_nl(4:end-2)*1e9,U_model_nl(4:end-2)-U_0_exp_nl, '--','LineWidth',3,'Color','r', 'DisplayName',  'Non -linear fitting')

 
box on
 

if aa==5
    set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
else
    set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25,'Yticklabel',[]);
end

xlabel('$x( \rm nm)$','Interpreter','Latex', 'FontSize',30)


if aa==5
    ylabel('$U(k_{\rm B} T)$','Interpreter','Latex','FontSize',30)
end
ylim([-0.5, 10])
 
hold off


xwi = 400;    % width of the plot square
bx1 = 140;     % extra space at the left
bx2 = 20;     % extra space at the right

Xpix = 3*xwi+bx1+3*bx2;;  % total

ywi = 300;    % length riquadro con funzione
by1 = 75;     % extra space below
by2 = 30;     % extra space up

Ypix = 2*by1+2*ywi+3*by2+20;  % larghezza figura in pixel
axes('Position',[(0) 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
hold on

xt = [bx1-110,bx1+xwi+bx2,bx1+2*(xwi+bx2),bx1-85,bx1+xwi+bx2,bx1+2*(xwi+bx2)];
yt = [ 2.6*by1+2*ywi+by2, 2.6*by1+2*ywi+by2, 2.6*by1+2*ywi+by2, ywi+3*by2+50, ywi+3*by2+50, ywi+3*by2+50];
str = {'\bf a','\bf b','\bf c', '\bf d','\bf e','\bf f'};
text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])


end
