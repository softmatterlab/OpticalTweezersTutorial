close all
clear all
load('../data/data_heating_rates.mat')

col3=[0.00 0.45 0.74; 216/255 82/255 24/255; 237/255 177/255 32/255; 128/255 128/255 128/255 ];

xwi = 390;
bx1 = 130;     % extra space at the left
bx2 = 100;     % extra space at the right
bxm = 20;     % extra space at the midle


Xpix=1400;   % total width 
ywi = 600;    % length frame with function

by1 = 100;     % extra space below
by2 = 30;     % extra space up
bym = 170      %extra space at the midle
Ypix = by1+ywi+2*by2;  % high in pixel


figure('Position',[100 -200 Xpix Ypix]);
axes('Position',[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix); 

% p(1)=errorbar(pressure,Gamma_ref, diff_Gamma_ref, diff_Gamma_ref, diff_pressure, diff_pressure, '.','Color', col3(1,:),'LineWidth',1 )
p(1)=scatter(pressure,Gamma_ref,100, 's','filled','MarkerEdgeColor', col3(1,:),'LineWidth',1 )
hold on

ft =fittype('x*a');
f1=fit(pressure', (Gamma_ref+diff_Gamma_ref/2)', ft, 'Weights', diff_Gamma_ref);

p(3)=plot(pressure, f1.a*pressure, '--','Color', col3(1,:),'LineWidth',  2)

% p(2)=errorbar(pressure,Gamma_co2, diff_Gamma_co2, diff_Gamma_co2, diff_pressure, diff_pressure, '.','Color', col3(2,:),'LineWidth',1)
p(2)=scatter(pressure,Gamma_co2,100, 's','filled','MarkerEdgeColor', col3(2,:),'LineWidth',1 )
f2=fit(pressure', (Gamma_co2+diff_Gamma_co2/2)', ft, 'Weights', diff_Gamma_co2);

p(4)=plot(pressure, f2.a*pressure,'--','Color', col3(2,:),'LineWidth', 2)
box on
LL=legend([p(1)  p(2) p(3) p(4)],'no heating', 'laser heating', 'fit', 'fit')
LL.Box = 'off';
LL.FontSize = 18;
LL.Position= [0.070 0.62 0.2164 0.2740]
LL.Interpreter='latex';

xlabel('Pressure $P_\mathrm{gas}$ (mbar)','Interpreter', 'latex')
ylabel('Heating Rate $\Gamma / (2\pi)$ (K/s)', 'Interpreter', 'latex')

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25, 'XMinorTick','on',...
      'yscale', 'log','xscale', 'log','Xtick',[1e-6 1e-3],'ylim', [1e-1, 3e3],'xlim',[1e-7, 1e-2])


%%
 load('../data/data_heating_rates_power_1e-05.mat')
axes('Position',[2*bx1+xwi 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix); 

errorbar(power*1000,Gamma_power, Gamma_power_error, Gamma_power_error, power*100, power*100, '.',...
    'Color', [0.99, 0.56, 0.2 ],'LineWidth',2)

xlabel('Laser Power $P$ (mW)','Interpreter', 'latex')
ylabel('Heating Rate $\Gamma / (2\pi)$ (K\ s$^{-1}$)','Interpreter', 'latex')

% txt = '$P_{\rm gas}=1x10^5$ (mbar)','Interpreter', 'Latex')
text(0,7.5,'$P_{\rm gas}$=1x10$^{-5}$ (mbar)','Interpreter','latex','FontSize',18)
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25, 'XMinorTick','on','xlim',[-1, 15])


%%
load('../data/data_frequency_shif.mat')
N=length(T);

colorMap = [linspace(1,0,N)', zeros(N,1),linspace(0,1,N)'];

axes('Position',[3*bx1+2*xwi 0 xwi/2 0]/Xpix + [0 by1 0 ywi]/Ypix); 


hold on 
% for i=1:N-1
% p(1)=scatter(T(i), delta_freq(i),100,'o','MarkerEdgeColor',colorMap(i,:),'LineWidth', 1)
% hold on
% end

p(1)=scatter(T, delta_freq,80,'o','MarkerEdgeColor','r','LineWidth', 0.5)
p(2)=plot(T, linear_fit, 'LineWidth', 2, 'color', 'k')

xlabel('$T_{\rm bulk}$ (K)','Interpreter', 'latex')
ylabel('$\Delta\Omega/\Omega_0$ (\%)', 'Interpreter', 'latex')
LL=legend([p(1) p(2)],'calc.','fit')
LL.Box = 'off';
LL.Position= [0.78 0.78 0.2164 0.1]
LL.FontSize = 18;
LL.Interpreter='latex';


box on
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25, 'XMinorTick','on','xlim',[300, 1600],'ylim',[0, 1.8])

%%
axes('Position',[0 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix); 


xt = [bx1,2*bx1+xwi,3*bx1+2*xwi];
yt = [by1+ywi+by2,by1+ywi+by2,by1+ywi+by2];
    
 str = {'\bf a','\bf b','\bf c'};


text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])


% saveas(gcf,'Fig46.eps','epsc')
% saveas(gcf,'Fig46.fig')