clear all
close all



xwi = 500;
bx1 = 130;     % extra space at the left
bx2 = 100;     % extra space at the right
bxm = 20;     % extra space at the midle


Xpix=1400;   % total width 
ywi = 600;    % length frame with function

by1 = 100;     % extra space below
by2 = 30;     % extra space up
bym = 170      %extra space at the midle
Ypix = by1+ywi+2*by2;  % high in pixel
% 
 load('Fig43_a.mat')
 figure('Position',[10 20 Xpix Ypix]);
 axes('Position',[bx1 0 xwi/2 0]/Xpix + [0 by1 0 ywi]/Ypix); 
  
 scatter (x1,y1)
 
  box on
 set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01],...
     'XMinorTick','on','ylim', [0.1,0.35],'xlim', [-80 -20],'Xtick',[-60 -40 ])
 %'xlim', [50,5000],'ylim', [1,300],'yscale', 'log', 'xscale', 'log'
 
axes('Position',[bx1+xwi/2+bxm 0 xwi/2 0]/Xpix + [0 by1 0 ywi]/Ypix); 

  scatter (x2,y2)
  
 box on
 
  set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01],...
     'XMinorTick','on', 'Yticklabels',{ },'ylim', [0.1,0.35],'xlim', [20,80], 'Xtick',[40 60 ])
 
%%
load('Fig43_b.mat')
axes('Position',[2*bx1+bxm+xwi 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix); 

scatter (x,y)
hold on
scatter (x,yy)
plot(xm,ym)
box on

 
 set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01],...
     'XMinorTick','on','xlim', [50,5000],'ylim', [1,300],'yscale', 'log', 'xscale', 'log')
 
%  
%    x1=xm(:,1);y1=xm(:,2);
%   x2=xm(:,3);y2=xm(:,4);
%   save('Fig43_a.mat','x1','y1','x2','y2')