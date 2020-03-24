clear all
close all
clc

load('../data/data_NoiseFloor.mat')

iname2 = '../data/TimeTrace_P=50mBar_fdr=135kHz_Vdr=30.2V_PSDdrFit.fig';

fig135 = openfig(iname2);
axesObjs135    = get(fig135, 'Children');  %axes handles
dataObjs135    = get(axesObjs135, 'Children'); %handles to low-level graphics objects in axes
xdata135    = dataObjs135(2).XData;
ydata135    = dataObjs135(2).YData;
xfit135     = dataObjs135(1).XData;
yfit135     = dataObjs135(1).YData;

col3=[0.00,0.45,0.74];


xwi = 1200;    % width of the plot square
bx1 = 170;     % extra space at the left
bx2 = 100;     % extra space at the right
bxm = 20;     % extra space at the midle


Xpix=1400;   % total width
ywi = 800;    % length frame with function

by1 = 100;     % extra space below
by2 = 30;     % extra space up
bym = 30      %extra space at the midle
Ypix = by1+ywi+by2;  % width in pixel


figure('Position',[10 20 Xpix Ypix]);
axes('Position',[bx1 0 xwi 0]/Xpix + [0 by1  0 ywi]/Ypix);

semilogy(xdata135,ydata135,'LineWidth', 2.5,'Color',col3)
hold on
semilogy(xfit135,yfit135,'--','LineWidth', 2.5,'Color','r')
plotNoise = semilogy(data_NoiseFloor(1,:)/1000,data_NoiseFloor(2,:),'LineWidth', 2.5,'Color',[0.5 0.5 0.5]);


xlabel('Frequency $\omega/ 2\pi$ (kHz)', 'Interpreter', 'Latex')
ylabel({'PSD','${S}_{v}(\omega)$ (bit$^2$ Hz$^{-1}$)'}, 'Interpreter', 'Latex')
text(140,40,'$S_v^{\rm th}(\omega_{\rm dr})$','Interpreter','latex','FontSize',30)
plot([100,150],[22.83,22.83],'--k', 'LineWidth',2.5)
text(140,700,'$S_v(\omega_{\rm dr})$','Interpreter','latex','FontSize',30)
plot([100,150],[1264,1264],'--k', 'LineWidth',2.5)
LL= legend ({'Data','Fit', 'Noise'},'Interpreter','latex','Box','off','Position',[0.17 0.69 0.1 0.2])
LL.FontSize = 18


% breakyaxis([0.02 50]);
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25, 'XMinorTick','on','xlim', [100,150],'ylim', [2e-3,3e3])%,...
   % 'Ytick',[1e-3 1e0 1e3 ]);


axes('Position',[bx1+300 0 0.4*xwi 0]/Xpix + [0 by1+160  0 0.4*ywi]/Ypix);

% plot(xdata135,ydata135,'-o','filled','MarkerEdgeColor',col3,'MarkerFaceColor',col3)
plot(xdata135,ydata135,'-o','MarkerEdgeColor',col3,'MarkerFaceColor',col3,'MarkerSize', 3)
box on
hold on
plot(xfit135,yfit135,'--','LineWidth', 2.5,'Color','r')
xlabel('Frequency $\omega/ 2\pi$ (kHz)', 'Interpreter', 'Latex')
ylabel({'PSD','${S}_{v}(\omega)$ (bit$^2$ Hz$^{-1}$)'}, 'Interpreter', 'Latex')



set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',18,'TickLength',[0.02, 0.01], 'XMinorTick','on','xlim', [134.8,135.2],'ylim', [8e0,2e3],...
    'Ytick',[1e1 1e2 1e3 ], 'yscale', 'log');



saveas(gcf,'Fig44.eps','epsc')
saveas(gcf,'Fig44.fig')
