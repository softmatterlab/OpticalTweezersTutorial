clear all
close all
set(0, 'defaultFigureRenderer', 'painters')

clear all

col3=[0.00,0.45,0.74];


xwi = 600;    % width of the plot square
xw2 = 300;
bx1 = 80;     % extra space at the left
bx2 = 100;     % extra space at the right
bxm = 30;     % extra space at the midle


Xpix=xwi+xw2+bx1+bxm+bxm;   % total width 
ywi = 400;    % length frame with function

by1 = 100;     % extra space below
by2 = 30;     % extra space up
bym = 30      %extra space at the midle
Ypix = by1+ywi+by2+bym;  % width in pixel


figure('Position',[10 20 Xpix Ypix]);

axes('Position',[bx1+xw2+bxm 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  

load('../data/data_clustering2.mat')
bar(x,y, 'FaceColor',col3 ,'EdgeColor','k','BarLayout','stacked')
xlabel('Effecttive Brightness B$_{\rm {eff}}$', 'Interpreter', 'Latex')
ylabel({'Counts'}, 'Interpreter', 'Latex')
text(200,45,'1 particle','Interpreter','latex','FontSize',25)
text(1000,24,'2 particle','Interpreter','latex','FontSize',25)
text(1950,18,'3 particle','Interpreter','latex','FontSize',25)
text(2800,10,'N$>$3','Interpreter','latex','FontSize',25)

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','on','xlim', [0,3500],'ylim', [0,50])
%  'Xtick',[1e-5  1e0 1e5 ],'Ytick',[1e2 1e4 1e6]);

% histogram(y,10,'DisplayStyle', 'bar', 'Binwidth',1.0,'Normalization', 'pdf','LineWidth',2, 'FaceColor',col3 ,'EdgeColor','k')

axes('Position',[bx1+xw2+bxm+310 0 0.45*xwi 0]/Xpix + [0 by1+245 0 ywi/3]/Ypix); 


load('../data/data_clustering.mat')

bar(x,y, 'FaceColor',col3 ,'EdgeColor','k','BarLayout','stacked')
xlabel('B$_{\rm {eff}}$', 'Interpreter', 'Latex')
ylabel({'Counts'}, 'Interpreter', 'Latex')

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',18,'TickLength',[0.02, 0.01], 'XMinorTick','on','xlim', [0,35000],'ylim', [0,100],...
    'Xtick',[0 1.5e4 3e4 ])

axes('Position',[0 0 xw2 0]/Xpix + [0 by1 0 ywi]/Ypix); 

img = imread('../figures/Fig40a.jpg');
% image('CData',img,'XData',[4 600],'YData',[0 500]);
image(img,'XData',[4 600],'YData',[0 500]);

set(gca, 'Xtick',[ ],'XMinorTick','off','YMinorTick','off','Ytick',[ ])

axes('Position',[0 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix); 


hold on

xt = [0,0.48*xw2,bx1+bxm+xw2];
yt = [by1+ywi+bym,0.7*ywi,by1+ywi+bym];
    
 str = {'\bf a','\bf b','\bf c'};


text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])


% saveas(gcf,'Fig40.eps','epsc')
% saveas(gcf,'Fig40.fig')