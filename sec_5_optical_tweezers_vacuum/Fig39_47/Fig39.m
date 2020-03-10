clear all
close all
DATA=load('data_pressure.mat')

DATAT=struct2table(DATA);


xwi = 1220;
bx1 = 140;     % extra space at the left
bx2 = 100;     % extra space at the right
bxm = 20;     % extra space at the midle


Xpix=1400;   % total width 
ywi = 600;    % length frame with function

by1 = 100;     % extra space below
by2 = 30;     % extra space up
bym = 170      %extra space at the midle
Ypix = by1+ywi+by2;  % high in pixel


figure('Position',[100 -100 Xpix Ypix]);
axes('Position',[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix); 

for j=1:size(DATAT,2)-2
     datai=DATAT.(j);
     p(j)=semilogy(datai.Time, datai.Pressure, 'DisplayName', datai.name,'LineWidth', 2)
     
     hold on
end
datai=DATAT.(6);
     p(j)=semilogy(datai.Time, datai.Pressure,'--', 'color','k','DisplayName', datai.name,'LineWidth', 2)
     
     datai=DATAT.(7);
     p(j)=semilogy(datai.Time, datai.Pressure, ':', 'color','k','DisplayName', datai.name,'LineWidth', 2)

ylabel('$P_{\rm gas}$ (mbar)','Interpreter', 'Latex')
xlabel('$t$ (h)','Interpreter', 'Latex')

LL=legend
LL.Box = 'off';
LL.FontSize = 18;
LL.NumColumns = 1;
LL.Position= [0.650 0.5703 0.2164 0.3740]
  

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25, 'XMinorTick','on',...
    'xlim',[-0.5, 11],'ylim', [1e-6, 3e-3],    'yscale', 'log')


saveas(gcf,'Fig39.eps','epsc')
saveas(gcf,'Fig39.fig')