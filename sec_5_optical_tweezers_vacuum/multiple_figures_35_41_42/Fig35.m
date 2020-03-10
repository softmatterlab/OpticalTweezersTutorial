clear all
col3=[0.00,0.45,0.74];

xwi = 500;    % width of the plot square
bx1 = 180;     % extra space at the left
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

LL= legend ([p(1) p(2) p(3) p(4)],'p=1000 mbar','p=60 mbar','p=2.5 mbar','fit','Interpreter','latex','Box','off','Position',[260 250 1 1], 'boxoff')
LL.FontSize = 18;
LL.Box = 'off';
hold off
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','on','xlim', [1e0,2e2],'ylim', [5e-8,1e-2],...
 'Xtick',[1e0 1e1 1e2 ],'Ytick',[1e-7 1e-5 1e-3 ]);

saveas(gcf,'Fig35.eps','epsc')
saveas(gcf,'Fig35.fig')