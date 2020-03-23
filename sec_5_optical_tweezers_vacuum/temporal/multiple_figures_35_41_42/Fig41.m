%%FIGURES FOR RUBEN
clear all
close all

%FIG1
clear all

col3=[0.00,0.45,0.74];


xwi = 210;    % width of the plot square
bx1 = 110;     % extra space at the left
bx2 = 100;     % extra space at the right
bxm = 33;     % extra space at the midle


Xpix=1400;   % total width 
ywi = 300;    % length frame with function

by1 = 100;     % extra space below
by2 = 30;     % extra space up
bym = 30      %extra space at the midle
Ypix = by1+ywi+by2+bym;  % width in pixel


figure('Position',[10 20 Xpix Ypix]);

axes('Position',[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  

load('data1_vacuum.mat')

loglog(pressure, gamma_tot/(2*pi), 'LineWidth', 2.5, 'Color',col3);
hold on
loglog(pressure, gamma_tot_2_bath/(2*pi),'--', 'LineWidth', 2.5,'Color','r')
xlabel('pressure (mbar)','Interpreter','latex')
ylabel('damping rate $\gamma/(2\pi)$ (Hz)','Interpreter','latex')

hold off
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','on','xlim', [9e-9,1e6],'ylim', [1e1,2e6],...
 'Xtick',[1e-5  1e0 1e5 ],'Ytick',[1e2 1e4 1e6]);

axes('Position',[2*bx1+xwi+bxm 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  
loglog(pressure, Gamma/kB, 'LineWidth', 2.5,'Color',col3);
hold on
loglog(pressure, Gamma_2_bath/kB,'--', 'LineWidth', 2.5,'Color','r')
xlabel('pressure (mbar)','Interpreter','latex')
ylabel('heating rate $\Gamma / k_{\rm B}$ (Ks$^{-1}$)','Interpreter','latex')

hold off
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','on','xlim', [1e-8,1e6],'ylim', [1e-1,5e9],...
 'Xtick',[1e-5  1e0 1e5 ],'Ytick',[1e0 1e4 1e8]);

axes('Position',[3*bx1+2*xwi+2*bxm 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  
loglog(pressure, T_FB, 'LineWidth', 2.5,'Color',col3);
hold on
loglog(pressure,  T_FB_2_bath,'--', 'LineWidth', 2.5,'Color','r')
xlabel('pressure (mbar)','Interpreter','latex')
ylabel('COM temperature (K)','Interpreter','latex')
hold off
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','on','xlim', [1e-8,1e6],'ylim', [1e-3,1e3],...
 'Xtick',[1e-5 1e0 1e5 ], 'Ytick',[1e-3 1e0 1e3 ]);

axes('Position',[4*bx1+3*xwi+3*bxm 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  
semilogx(pressure, T_int, 'LineWidth', 2.5,'Color',col3);
hold on
xlabel('pressure (mbar)','Interpreter','latex')
ylabel('internal temperature (K)','Interpreter','latex')
hold off
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','on','xlim', [1e-8,1e6],'ylim', [200,1200],...
 'Xtick',[1e-5 1e0 1e5 ]);

axes('Position',[0 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix); 
hold on

xt = [bx1,2*bx1+xwi+bxm,3*bx1+2*xwi+2*bxm,4*bx1+3*xwi+3*bxm];
yt = [by1+ywi+by2,by1+ywi+by2,by1+ywi+by2,by1+ywi+by2];
    
 str = {'\bf a','\bf b','\bf c','\bf d'};


text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])

saveas(gcf,'Fig41.eps','epsc')
saveas(gcf,'Fig41.fig')