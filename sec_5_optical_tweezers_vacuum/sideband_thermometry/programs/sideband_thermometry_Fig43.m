clear all
close all

col3=[0.00,0.45,0.74];
col4=[216/255 82/255 24/255];

xwi = 550;
bx1 = 130;     % extra space at the left
bx2 = 100;     % extra space at the right
bxm = 20;     % extra space at the midle


Xpix=1400;   % total width 
ywi = 500;    % length frame with function

by1 = 100;     % extra space below
by2 = 30;     % extra space up
bym = 170      %extra space at the midle
Ypix = by1+ywi+2*by2;  % high in pixel
% 
 load('../data/Fig43_a.mat')
 delta=0.05
 id1 = find(x1> -60-delta & x1< -60+delta);
 id2 = find(x1> -40-delta & x1< -40+delta);
  xa = linspace(-80,-20,100);
 ya= linspace(0.13,0.13,100);
 figure('Position',[10 20 Xpix Ypix]);
 axes('Position',[bx1 0 xwi/2 0]/Xpix + [0 by1 0 ywi]/Ypix); 
 
 rectangle('Position',[-80,0.1,60,0.03],'FaceColor',[153 153 153]/255)
%  rectangle('Position',[-80,0.1,60,0.03],'FaceColor',[0 .5 .5],'EdgeColor','b',  'LineWidth',3)
  hold on
 plot (x1,y1,'LineWidth', 1.5, 'Color', 'k')

 plot([-60 -60],[0 0.5],'--','LineWidth', 1.5, 'color','k')
 plot([-40 -40],[0 0.5],'--','LineWidth', 1.5, 'color','k')

  plot([-80 -20],[0.13 0.13],'-','LineWidth', 1.5, 'color','k')

area(x1(id1:id2), y1(id1:id2),0.13,'LineStyle','none','faceColor', 'r','FaceAlpha' , 0.5)

 
 xlabel('$\Delta f$ (kHz)','Interpreter', 'Latex')
 ylabel('PSD (pm$^2$  Hz$^{-1}$)','Interpreter', 'Latex')

 text(-45,0.3,'$\propto \bar n + 1$','Interpreter','latex','FontSize',22,'BackgroundColor', 'w')
 
  box on
 set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01],...
     'XMinorTick','on','ylim', [0.1,0.35],'xlim', [-80 -20],'Xtick',[-60 -40 ])

 %%
 delta=0.05
  id1 = find(x2> 40-delta & x2< 40+delta);
 id2 = find(x2> 60-delta & x2< 60+delta);
 
 xa = linspace(20,80,100);

 
axes('Position',[bx1+xwi/2+bxm 0 xwi/2 0]/Xpix + [0 by1 0 ywi]/Ypix); 

rectangle('Position',[20,0.1,60,0.03],'FaceColor',[153 153 153]/255)
  hold on
  plot (x2,y2,'LineWidth', 1.5, 'Color','k')

  plot([60 60],[0 0.5],'--','LineWidth', 1.5, 'color','k')
  plot([40 40],[0 0.5],'--','LineWidth', 1.5, 'color','k')
  
  plot([20 80],[0.13 0.13],'-','LineWidth', 1.5, 'color','k')
  
  area(x2(id1:id2), y2(id1:id2),0.13,'LineStyle','none','faceColor', col3,'FaceAlpha' , 0.6)
  
  text(55,0.3,'$\propto \bar n$','Interpreter','latex','FontSize',22,'BackgroundColor', 'w')
  
 xlabel('$\Delta f$ (kHz)','Interpreter', 'Latex')
  
 box on
 
  set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01],...
     'XMinorTick','on', 'Yticklabels',{ },'ylim', [0.1,0.35],'xlim', [20,80], 'Xtick',[40 60 ])
 
%%
load('../data/Fig43_b.mat')
axes('Position',[2*bx1+bxm+xwi 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix); 
p(3)=plot(xm,ym,'--','Color', 'k','LineWidth', 2)
hold on
p(1)=scatter (x,y,75,'o','filled','MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth', 1)

p(2)=scatter (x,yy,75,'d','filled','MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth', 1)

box on

 xlabel('$\gamma_{\rm fb} / (2\pi)$ (kHz)','Interpreter', 'Latex')
 ylabel('$\bar n + 1/2 $','Interpreter', 'Latex')

LL=legend([p(1) p(2) p(3) ],'asymmetry', 'left sideband', 'model')
LL.Box = 'off';
LL.FontSize = 18;
LL.Position= [0.75 0.65 0.2164 0.2]
LL.Interpreter='latex';
 
 set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01],...
     'XMinorTick','on','xlim', [7e0,2e4],'ylim', [2e0,3e3],'yscale', 'log', 'xscale', 'log')
 %%
axes('Position',[0 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix); 


xt = [bx1,bx1+bxm+xwi/2,2*bx1+xwi+bxm];
yt = [by1+ywi+by2,by1+ywi+by2,by1+ywi+by2];
    
 str = {'\bf a','\bf b','\bf c'};


text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])


% saveas(gcf,'Fig42.eps','epsc')
% saveas(gcf,'Fig42.fig')