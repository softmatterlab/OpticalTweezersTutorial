%%FIGURES FOR RUBEN
clear all
close all

%FIG1
clear all

col3=[0.00,0.45,0.74];


xwi = 215;    % width of the plot square
bx1 = 120;     % extra space at the left
bx2 = 100;     % extra space at the right
bxm = 20;     % extra space at the midle


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

axes('Position',[2*bx1+xwi 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  
loglog(pressure, Gamma/kB, 'LineWidth', 2.5,'Color',col3);
hold on
loglog(pressure, Gamma_2_bath/kB,'--', 'LineWidth', 2.5,'Color','r')
xlabel('pressure (mbar)','Interpreter','latex')
ylabel('heating rate $\Gamma / k_{\rm B}$ (Ks$^{-1}$)','Interpreter','latex')

hold off
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','on','xlim', [1e-8,1e6],'ylim', [1e-1,5e9],...
 'Xtick',[1e-5  1e0 1e5 ],'Ytick',[1e0 1e4 1e8]);

axes('Position',[3*bx1+2*xwi+bxm 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  
loglog(pressure, T_FB, 'LineWidth', 2.5,'Color',col3);
hold on
loglog(pressure,  T_FB_2_bath,'--', 'LineWidth', 2.5,'Color','r')
xlabel('pressure (mbar)','Interpreter','latex')
ylabel('COM temperature (K)','Interpreter','latex')
hold off
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','on','xlim', [1e-8,1e6],'ylim', [1e-3,1e3],...
 'Xtick',[1e-5 1e0 1e5 ], 'Ytick',[1e-3 1e0 1e3 ]);

axes('Position',[4*bx1+3*xwi+2*bxm 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  
semilogx(pressure, T_int, 'LineWidth', 2.5,'Color',col3);
hold on
xlabel('pressure (mbar)','Interpreter','latex')
ylabel('internal temperature (K)','Interpreter','latex')
hold off
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','on','xlim', [1e-8,1e6],'ylim', [200,1200],...
 'Xtick',[1e-5 1e0 1e5 ]);

axes('Position',[0 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix); 
hold on

xt = [bx1,2*bx1+xwi,3*bx1+2*xwi+bxm,4*bx1+3*xwi+2*bxm];
yt = [by1+ywi+by2,by1+ywi+by2,by1+ywi+by2,by1+ywi+by2];
    
 str = {'\bf a','\bf b','\bf c','\bf d'};


text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])

saveas(gcf,'Fig41.eps','epsc')
saveas(gcf,'Fig41.fig')


%%
clear all
col3=[0.00,0.45,0.74];

xwi = 480;    % width of the plot square
bx1 = 200;     % extra space at the left
bx2 = 100;     % extra space at the right
bxm = 20;     % extra space at the midle


Xpix=700;   % total width 
ywi = 300;    % length frame with function

by1 = 100;     % extra space below
by2 = 20;     % extra space up
bym = 30      %extra space at the midle
Ypix = by1+ywi+by2;  % width in pixel

load('data2_vacuum_0.mat')

figure('Position',[10 20 Xpix Ypix]);

axes('Position',[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  

p(1)=loglog(frequency/1000, PSD/g0, 'LineWidth', 2.5, 'color',col3)
hold on 
p(4)=loglog(frequency/1000, Lorentz,'--', 'LineWidth', 2.5, 'Color', 'black')
plot([x0 x0],[ymin ymax-ymax/2],'--', 'LineWidth', 2.5, 'color',col3)

load('data2_vacuum_1.mat')
p(2)=loglog(frequency/1000, PSD/g0, 'LineWidth', 2.5, 'Color', ' [0.9000 0.400 0.1000]' )
hold on 
loglog(frequency/1000, Lorentz,'--', 'LineWidth', 2.5, 'Color', 'black')
plot([x0 x0],[ymin ymax-ymax/5],'--','LineWidth', 2.5, 'Color', ' [0.9000 0.400 0.1000]')

load('data2_vacuum_2.mat')
p(3)=loglog(frequency/1000, PSD/g0, 'LineWidth', 2.5, 'Color', '[0 0.7 0]')
hold on 
loglog(frequency/1000, Lorentz,'--', 'LineWidth', 2.5, 'Color', 'black')
plot([x0 x0],[ymin ymax],'--','LineWidth', 2.5, 'Color', 'r')
xlabel('Frequency $f$ (kHz)', 'Interpreter', 'Latex')
ylabel({'Normalized PSD';'$\hat{S}_{vv}(f)/g$ (bit$^2$ {Hz}$^{-2})$'}, 'Interpreter', 'Latex')

LL= legend ([p(1) p(2) p(3) p(4)],'1000 mbar','60 mbar','2.5 mbar','fit','Interpreter','latex','Box','off','Position',[260 250 1 1], 'boxoff')
LL.FontSize = 18;
LL.Box = 'off';
hold off
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','on','xlim', [1e0,2e2],'ylim', [5e-8,1e-2],...
 'Xtick',[1e0 1e1 1e2 ],'Ytick',[1e-7 1e-5 1e-3 ]);

saveas(gcf,'Fig35.eps','epsc')
saveas(gcf,'Fig35.fig')


%%
clear all
close all
col3=[0.00,0.45,0.74];

xwi = 500;    % width of the plot square
bx1 = 180;     % extra space at the left
bx2 = 60;     % extra space at the right
bxm = 20;     % extra space at the midle


Xpix=700;   % total width 
ywi = 300;    % length frame with function

by1 = 100;     % extra space below
by2 = 20;     % extra space up
bym = 20      %extra space at the midle
Ypix = by1+2*ywi+bym+by2;  % width in pixel


load('data3_vacuum.mat')

figure('Position',[500 20 Xpix Ypix]);

axes('Position',[bx1 0 xwi 0]/Xpix + [0 by1+ywi+bym  0 ywi]/Ypix);  
 
% p(1)=semilogy(f1, PSD1, 'LineWidth', 2.5, 'color',col3)

p(1)=scatter(f1, PSD1,25,'o','filled','MarkerEdgeColor',col3,'MarkerFaceColor',col3)


hold on
p(2)=plot([xx0 xx0],[1e0 1e3],'--','LineWidth', 1.5, 'color',col3)
plot([xx1 xx1],[1e0 1e3],'--','LineWidth', 1.5, 'color',col3)
plot([xx2 xx2],[1e0 1e3],'--','LineWidth', 1.5, 'color',col3)
ylabel({'PSD','$\hat{S}_{qq}$ (nm$^2$ Hz$^{-2}$)'}, 'Interpreter', 'Latex')

%{'PSD';'$\hat{S}_{qq}(f)/g \ (bit^2\rm\ {Hz}^{-2})$'}
box on
area(f1, PSD1,'LineStyle','none','faceColor', col3,'FaceAlpha' , 0.4)
% a.FaceAlpha = 0.5 

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','on','xlim', [20,60],'ylim', [4e0,1e3],...
  'Xtick',[20 30 40 50 60 ],'Xticklabel',[],'Ytick',[1e1 1e2],'yscale', 'log');


axes('Position',[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  


scatter(f2, PSD2,25,'o','filled','MarkerEdgeColor',col3,'MarkerFaceColor',col3)

hold on
plot([XV25000 XV25000],[1e-2 1e3],'--','LineWidth', 1.5, 'color',col3)
plot([XV45000 XV45000],[1e-2 1e3],'--','LineWidth', 1.5, 'color',col3)
plot([XV50000 XV50000],[1e-2 1e3],'--','LineWidth', 1.5, 'color',col3)
ylabel({'PSD','$\hat{S}_{vv}$ (bit$^2$ Hz$^{-2}$)'}, 'Interpreter', 'Latex')


LL= legend ([p(1) p(2)],'z detector signal','drive', 'Interpreter', 'Latex','Box','off','Position',[275 500 1 1], 'boxoff')
LL.FontSize = 18;
LL.Box = 'off';

xlabel('frequency $f$ (kHz)',  'Interpreter', 'Latex')
box on
area(f2, PSD2,'LineStyle','none','faceColor', col3,'FaceAlpha' , 0.4)
 set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','on','xlim', [20,60],'ylim', [2e-3,3e2],'yscale', 'log',...
   'Xtick',[20 40 60],'Ytick',[1e-2 1e0 1e2 ]);
axes('Position',[0 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix); 


hold on

xt = [0,0];
yt = [by1+2*ywi+bym,by1+ywi];
    
 str = {'\bf a','\bf b'};


text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])


saveas(gcf,'Fig42.eps','epsc')
saveas(gcf,'Fig42.fig')
