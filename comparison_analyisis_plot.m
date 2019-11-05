
close all
clear all
load('results_averaged_Data_positions_Fig9_1P2_S.mat20191104T034638.mat')

xwi = 500;    % width of the plot square
bx1 = 110;     % extra space at the left
bx2 = 40;     % extra space at the right

Xpix = 1*xwi+1*bx1+1*bx2;  % total

ywi = 250;    % length riquadro con funzione
by1 = 100;     % extra between plots
by2 = 150;     % extra below
Ypix = 4*by1+4*ywi+by2;  % larghezza figura in pixel
%number of bins of the histogram, if not set default is 50  1
cyan=[0/255,200/255,158/255];
yellow=[233/255, 189/255, 21/255];
gray_blue=[0.1 0.5 0.7];
red_wine=[0.5 0.05 0.3];
figure('Position',[10 20 Xpix Ypix]);
positioninthefig1=[bx1 0 xwi 0]/Xpix + [0 3*by1+3*ywi+by2 0 ywi]/Ypix ;

axes( 'Position',positioninthefig1);  % fa in modo di centrare il riquadro degli assi nella posizione voluta

%load('results_comparison.mat')
semilogx(dnn, 1e6*kk_pot,'DisplayName', 'Potential', 'Color', 'red', 'LineWidth', 2)
hold on
semilogx(dnn, 1e6*kk_eq, '--','DisplayName', 'Equipartition', 'Color', 'blue', 'LineWidth', 2)
semilogx(dnn, 1e6*kk_acf, 'DisplayName', 'ACF', 'Color', cyan, 'LineWidth',2)
semilogx(dnn, 1e6*kk_psd, 'DisplayName', 'PSD', 'Color', yellow, 'LineWidth', 2)
semilogx(dnn, 1e6*kk_msd, 'DisplayName', 'MSD', 'Color', red_wine, 'LineWidth', 2)
semilogx(dnn, 1e6*kk_bay,'DisplayName', 'BAYESIAN', 'Color', gray_blue, 'LineWidth', 2)
semilogx(dnn, 1e6*kk_forma, '--','DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 2)



legend('FontSize', 15, 'Interpreter', 'Latex', 'NumColumns',2)
legend boxoff 
%xlabel('$N_s$', 'Interpreter', 'Latex', 'FontSize',30)
ylab=ylabel('$\kappa( p \rm N  \mu \rm m^{-1})$', 'Interpreter', 'Latex', 'FontSize',25)
set(ylab, 'Units', 'Normalized', 'Position', [-0.08, 0.5, 0]);
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
%set(gca,'XTickLabel',)
ylim([12 21])
xlim([100 330000])
hold off
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25);

%%
positioninthefig2=[bx1 0 xwi 0]/Xpix + [0 2*by1+2*ywi+by2 0 ywi]/Ypix ;

axes( 'Position',positioninthefig2);  % fa in modo di centrare il riquadro degli assi nella posizione voluta

semilogx(dnn,100* sigma_kk_pot./kk_pot,'DisplayName', 'Potential', 'Color', 'red', 'LineWidth', 2)
hold on
semilogx(dnn, 100*sigma_kk_eq./kk_eq,'--', 'DisplayName', 'Equipartition', 'Color', 'blue', 'LineWidth', 2)
semilogx(dnn, 100*sigma_kk_acf./kk_acf, 'DisplayName', 'ACF', 'Color', cyan, 'LineWidth', 2)
semilogx(dnn, 100*sigma_kk_psd./kk_psd, 'DisplayName', 'PSD', 'Color', yellow, 'LineWidth', 2)
semilogx(dnn, sigma_kk_msd./kk_msd, 'DisplayName', 'MSD', 'Color', red_wine, 'LineWidth', 2)

semilogx(dnn, 100*sigma_kk_bay./kk_bay,'DisplayName', 'BAYESIAN',  'Color', gray_blue, 'LineWidth', 2)
semilogx(dnn, 100*sigma_kk_forma./kk_forma, '--','DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 2)



ylab=ylabel('$\Delta \kappa/\kappa(\%)$', 'Interpreter', 'Latex', 'FontSize',25)
set(ylab, 'Units', 'Normalized', 'Position', [-0.08, 0.5, 0]);

set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);

ylim([-0.1 3])
xlim([100 330000])
hold off
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25);


%%
positioninthefig3=[bx1 0 xwi 0]/Xpix + [0 1*by1+1*ywi+by2 0 ywi]/Ypix ;

axes( 'Position',positioninthefig3); 


semilogx(dnn, DD_acf*1e12, 'DisplayName', 'ACF', 'Color', cyan, 'LineWidth', 2)
hold on
semilogx(dnn, DD_psd*1e12, 'DisplayName', 'PSD', 'Color', yellow, 'LineWidth', 2)
semilogx(dnn, DD_msd*1e12, 'DisplayName', 'MSD', 'Color', red_wine, 'LineWidth', 2)
semilogx(dnn, DD_bay*1e12, 'DisplayName', 'BAYESIAN', 'Color', gray_blue, 'LineWidth', 2)
semilogx(dnn, DD_forma*1e12, '--','DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 2)



semilogx(dnn, 0.208627*ones(size(dnn)), '--', 'DisplayName', 'Theoretical value', 'Color','black' , 'LineWidth', 2)
legend('FontSize', 15, 'Interpreter', 'Latex', 'NumColumns',2)
legend boxoff 
ylab=ylabel('$D(\mu \rm m^2/s)$', 'Interpreter', 'Latex', 'FontSize',25)
set(ylab, 'Units', 'Normalized', 'Position', [-0.08, 0.5, 0]);
ylim([0.17 0.35])
xlim([100 330000])

yticks([0, 0.1 0.2, 0.3])
box on 
hold off
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25);



%%
positioninthefig4=[bx1 0 xwi 0]/Xpix + [0 by2 0 ywi]/Ypix ;

axes( 'Position',positioninthefig4);  % fa in modo di centrare il riquadro degli assi nella posizione voluta


semilogx(dnn,100*sigma_DD_acf./ DD_acf, 'DisplayName', 'ACF', 'Color', cyan, 'LineWidth', 2)
hold on
semilogx(dnn, 100*sigma_DD_psd./DD_psd, 'DisplayName', 'PSD', 'Color', yellow, 'LineWidth', 2)

semilogx(dnn, 100*sigma_DD_msd./DD_msd, 'DisplayName', 'MSD', 'Color', red_wine, 'LineWidth', 2)
semilogx(dnn,100*sigma_DD_bay./DD_bay, 'DisplayName', 'BAYESIAN', 'Color', gray_blue, 'LineWidth', 2)
semilogx(dnn, 100*sigma_DD_forma./DD_forma, '--','DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 2)



xlabel('$N_s$', 'Interpreter', 'Latex', 'FontSize',25)
ylab=ylabel('$\Delta D/D(\%)$', 'Interpreter', 'Latex', 'FontSize',23)
set(ylab, 'Units', 'Normalized', 'Position', [-0.08, 0.5, 0]);
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
ylim([-0.1 2])

xlim([100 330000])
hold off
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25);


axes('Position',[(10) 0 Xpix 0]/Xpix + [0 20 0 Ypix]/Ypix);
hold on

yt = [4*by1+4*ywi+by2/2,3*by1+3*ywi+by2/2, 2*by1+2*ywi+by2/2,by1+ywi+by2/2]/Ypix;
xt = ones(size(yt))*(bx1/20)/Xpix;
str = {'\bf a','\bf b','\bf c', '\bf d'};

text(xt,yt,str,'Interpreter','Latex','FontSize',34)


hold off
axis off
