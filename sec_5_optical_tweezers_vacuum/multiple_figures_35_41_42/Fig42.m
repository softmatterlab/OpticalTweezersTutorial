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

figure('Position',[500 -200 Xpix Ypix]);

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
