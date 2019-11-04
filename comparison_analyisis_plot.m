

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
figure('Position',[10 20 Xpix Ypix]);
positioninthefig1=[bx1 0 xwi 0]/Xpix + [0 3*by1+3*ywi+by2 0 ywi]/Ypix ;

axes( 'Position',positioninthefig1);  % fa in modo di centrare il riquadro degli assi nella posizione voluta

%load('results_comparison.mat')
semilogx(dnn, 1e6*kk_pot, 'DisplayName', 'Potential', 'Color', 'red', 'LineWidth', 3)
hold on
semilogx(dnn, 1e6*kk_eq, 'DisplayName', 'Equipartition', 'Color', 'blue', 'LineWidth', 3)
semilogx(dnn, 1e6*kk_acf, 'DisplayName', 'ACF', 'Color', 'cyan', 'LineWidth', 3)
semilogx(dnn, 1e6*kk_psd, 'DisplayName', 'PSD', 'Color', 'yellow', 'LineWidth', 3)
semilogx(dnn, 1e6*kk_msd, 'DisplayName', 'MSD', 'Color', 'black', 'LineWidth', 3)

semilogx(dnn, 1e6*kk_forma, 'DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 3)
hold on
semilogx(dnn, 1e6*kk_bay, '--','DisplayName', 'BAYESIAN', 'Color', [0.1 0.5 0.7], 'LineWidth', 3)
%
legend('FontSize', 10)
%xlabel('$N_s$', 'Interpreter', 'Latex', 'FontSize',30)
ylabel('$k( p \rm N  \mu \rm m^{-1})$', 'Interpreter', 'Latex', 'FontSize',25)
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
%set(gca,'XTickLabel',)
ylim([8 20])
xlim([0 330000])
hold off
%%
positioninthefig2=[bx1 0 xwi 0]/Xpix + [0 2*by1+2*ywi+by2 0 ywi]/Ypix ;

axes( 'Position',positioninthefig2);  % fa in modo di centrare il riquadro degli assi nella posizione voluta

semilogx(dnn,100* sigma_kk_pot./kk_pot, 'DisplayName', 'Potential', 'Color', 'red', 'LineWidth', 3)
hold on
semilogx(dnn, 100*sigma_kk_eq./kk_eq, 'DisplayName', 'Equipartition', 'Color', 'blue', 'LineWidth', 3)
semilogx(dnn, 100*sigma_kk_acf./kk_acf, 'DisplayName', 'ACF', 'Color', 'cyan', 'LineWidth', 3)
semilogx(dnn, 100*sigma_kk_psd./kk_psd, 'DisplayName', 'PSD', 'Color', 'yellow', 'LineWidth', 3)
semilogx(dnn, sigma_kk_msd./kk_msd, 'DisplayName', 'MSD', 'Color', 'black', 'LineWidth', 3)


semilogx(dnn, 100*sigma_kk_forma./kk_forma, 'DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 3)

semilogx(dnn, 100*sigma_kk_bay./kk_bay,'--' ,'DisplayName', 'BAYESIAN',  'Color', [0.1 0.5 0.7], 'LineWidth', 3)

%xlabel('$N_s$', 'Interpreter', 'Latex', 'FontSize',30)
ylabel('$\Delta k/k(\%)$', 'Interpreter', 'Latex', 'FontSize',25)
legend('FontSize', 10)
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
%set(gca,'XTickLabel',[])
ylim([-0.1 5])
xlim([0 330000])
hold off

%%
positioninthefig3=[bx1 0 xwi 0]/Xpix + [0 1*by1+1*ywi+by2 0 ywi]/Ypix ;

axes( 'Position',positioninthefig3);  % fa in modo di centrare il riquadro degli assi nella posizione voluta


semilogx(dnn, DD_acf*1e12, 'DisplayName', 'ACF', 'Color', 'cyan', 'LineWidth', 3)
hold on
semilogx(dnn, DD_psd*1e12, 'DisplayName', 'PSD', 'Color', 'yellow', 'LineWidth', 3)
semilogx(dnn, DD_msd*1e12, 'DisplayName', 'MSD', 'Color', 'black', 'LineWidth', 3)

semilogx(dnn, DD_forma*1e12, 'DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 3)

semilogx(dnn, DD_bay*1e12,'--', 'DisplayName', 'BAYESIAN', 'Color', [0.1 0.5 0.7], 'LineWidth', 3)

semilogx(dnn, 0.208627*ones(size(dnn)), '--', 'DisplayName', 'Theoretical value', 'Color', [0.5 0.05 0.3])
legend('FontSize', 10)

%xlabel('$N_s$', 'Interpreter', 'Latex', 'FontSize',30)
ylabel('$D(\mu m^2/s)$', 'Interpreter', 'Latex', 'FontSize',25)
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
%set(gca,'XTickLabel',[])
ylim([0.1 0.3])
xlim([0 330000])
hold off

%%
positioninthefig4=[bx1 0 xwi 0]/Xpix + [0 by2 0 ywi]/Ypix ;

axes( 'Position',positioninthefig4);  % fa in modo di centrare il riquadro degli assi nella posizione voluta


%semilogx(dnn, D_eq*1e6, 'DisplayName', 'Equipartition', 'Color', 'blue', 'LineWidth', 3)
semilogx(dnn,100*sigma_DD_acf./ DD_acf, 'DisplayName', 'ACF', 'Color', 'cyan', 'LineWidth', 3)
hold on
semilogx(dnn, 100*sigma_DD_psd./DD_psd, 'DisplayName', 'PSD', 'Color', 'yellow', 'LineWidth', 3)
semilogx(dnn, 100*sigma_DD_msd./DD_msd, 'DisplayName', 'MSD', 'Color', 'black', 'LineWidth', 3)
semilogx(dnn, 100*sigma_DD_forma./DD_forma, 'DisplayName', 'FORMA', 'Color', 'magenta', 'LineWidth', 3)

semilogx(dnn,100*sigma_DD_bay./DD_bay,'--', 'DisplayName', 'BAYESIAN', 'Color', [0.1 0.5 0.7], 'LineWidth', 3)
legend('FontSize', 10)
xlabel('$N_s$', 'Interpreter', 'Latex', 'FontSize',25)
ylabel('$\Delta D/D(\%)$', 'Interpreter', 'Latex', 'FontSize',23)
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
ylim([-0.1 4])

xlim([0 330000])
hold off