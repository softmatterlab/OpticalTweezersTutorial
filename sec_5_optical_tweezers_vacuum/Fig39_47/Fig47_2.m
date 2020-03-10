clear all
clc
close all
set(0, 'defaultFigureRenderer', 'painters')
DATA=load('data_feeding_cooling_vacuum.mat');
DATAT=struct2table(DATA);


col3=[0.00 0.45 0.74;, 216/255 82/255 24/255; 237/255 177/255 32/255; 128/255 128/255 128/255 ];


xwi = 600;    % width of the plot square
xw2 = 500;
bx1 = 100;     % extra space at the left
bx2 = 100;     % extra space at the right
bxm = 150;     % extra space at the midle


Xpix=1400;   % total width 
ywi = 400;    % length frame with function

by1 = 100;     % extra space below
by2 = 30;     % extra space up
bym = 30      %extra space at the midle
Ypix = by1+ywi+by2+bym;  % width in pixel


figure('Position',[10 20 Xpix Ypix]);

axes('Position',[bx1 0 xw2 0]/Xpix + [0 by1 0 ywi]/Ypix);  

for j=1:size(DATAT,2)
    datai=DATAT.(j);
     Ft=lorentzfit(datai.freq,datai.psd);
    if j <size(DATAT,2)
 f(j)=plot(datai.freq, Ft,  'LineWidth', 2, 'color', 'k');
    end  
 hold on
    p(j)=scatter(datai.freq, datai.psd,55,'o','MarkerEdgeColor',[col3(j,:)],'LineWidth', 1);
    
    
end

%  text(142.4,3.4,'$\gamma_{\rm fb} =$','Interpreter','latex','FontSize',22)


ylabel('PSD $\tilde{S}^{\rm ol}_{\rm yy}$ (pm$^2$ Hz$^{-1}$)','Interpreter', 'Latex')
xlabel('$f$ (kHz)','Interpreter', 'Latex')
 
box on

 
 LL= legend ([p(1)  p(3) p(2) p(4) f(1)],'\gamma_{\rm fb} = 0.3 kHz','\gamma_{\rm fb} = 0.9 kHz',...
     '\gamma_{\rm fb} = 3.0 kHz','shotnoise','Lorentzian fit','Interpreter','latex','Box','off','Position',[155 280 1 1], 'boxoff')

% 
%  LL= legend ([p(1)  p(3) p(2) p(4) f(1)],'0.3 kHz','0.9 kHz',...
%      '3.0 kHz','shotnoise','fit','Interpreter','latex','Box','off','Position',[135 280 1 1], 'boxoff')
 
LL.FontSize = 18;
LL.Box = 'off';
  
 
 set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01],...
     'XMinorTick','on','xlim', [141,151],'ylim', [5e-1,4],'yscale', 'log')
 



 %%
 
 DATA=load('data_temperature.mat')
DATAT=struct2table(DATA);
axes('Position',[bx1+xw2+bxm 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix); 

datai=DATAT.(1);
p(1)=scatter(datai.gamma2pi, datai.T,85,'o','filled',    'MarkerFaceColor', [col3(2,:)]);
    hold on
    
datai=DATAT.(3);
p(2)=loglog(datai.gamma2pi, datai.T,'--', 'DisplayName', datai.name, 'LineWidth', 2, 'color','k');
    
    datai=DATAT.(4);
p(3)=scatter(datai.gamma2pi, datai.T,85,'o','filled',    'MarkerFaceColor', [col3(1,:)]);
 

datai=DATAT.(2);
loglog(datai.gamma2pi, datai.T, '--', 'DisplayName', datai.name, 'LineWidth', 2, 'color',col3(1,:));


datai=DATAT.(5); 
loglog(datai.gamma2pi, datai.T,'--',  'DisplayName', datai.name, 'LineWidth', 2, 'color','r');

datai=DATAT.(6);
loglog(datai.gamma2pi, datai.T,'--', 'DisplayName', datai.name, 'LineWidth', 2, 'color','k');


 LL= legend ([p(1) p(3) p(2) ],'p = 1.4 x 10^{-8} mbar','p = 1.2 x 10^{-7} mbar',...
     'model','Interpreter','latex','Box','off','Position',[620 120 1 1], 'boxoff')

% LL= legend ([p(1) p(3) p(2) ],'1.4x10^{-8} mbar','1.2x10^{-7} mbar',...
%      'model','Interpreter','latex','Box','off','Position',[525 120 1 1], 'boxoff')
LL.FontSize = 18;
LL.Box = 'off';

hold off
ylabel('$T_y$ (mK)','Interpreter', 'Latex')
xlabel('$\gamma_{\rm fb}/2\pi$ (kHz)','Interpreter', 'Latex')
box on
 set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'Xtick',[1e-2 1e-1 1e0 1e1],'Ytick',[1e-1 1e0 1e1],...
      'XMinorTick','on','xlim', [3e-3,1e1],'ylim', [5e-2,1e2],'yscale', 'log','xscale', 'log')



 %%
 



DATA=load('data_heterodyne_spectra_vacuum.mat')
DATAT=struct2table(DATA);
axes('Position',[bx1+xw2+bxm+230 0 0.4*xwi 0]/Xpix + [0 by1+250 0 0.35*ywi]/Ypix); 

% for j=1:size(DATAT,2)
%      datai=DATAT.(j);
%      Ft=lorentzfit(datai.freq,datai.psd);
%     if j <size(DATAT,2)
%  f(j)=plot(datai.freq, Ft,  'LineWidth', 2, 'color', 'k');
%     end  
%  hold on
%     p(j)=scatter(datai.freq, datai.psd,55,'o','MarkerEdgeColor',[col3(j,:)]);
%       
% end
% figure

 datai=DATAT.(1);
     Ft=lorentzfit(datai.freq,datai.psd);
   ms=35;
 f(1)=plot(datai.freq, Ft,  'LineWidth', 2, 'color', 'k');
   
 hold on
    p(1)=scatter(datai.freq, datai.psd,ms,'o','MarkerEdgeColor',[col3(1,:)]);
    
    
    datai=DATAT.(2);
    
    p(2)=scatter(datai.freq, datai.psd,ms,'o','MarkerEdgeColor',[col3(4,:)]);

datai=DATAT.(3);
     Ft=lorentzfit(datai.freq,datai.psd);
   
 f(3)=plot(datai.freq, Ft,  'LineWidth', 2, 'color', 'k');
   
 hold on
    p(3)=scatter(datai.freq, datai.psd,ms,'o','MarkerEdgeColor',[col3(2,:)]);
    
    datai=DATAT.(4);
     Ft=lorentzfit(datai.freq,datai.psd);
   
 f(4)=plot(datai.freq, Ft,  'LineWidth', 2, 'color', 'k');
   
 hold on
    p(4)=scatter(datai.freq, datai.psd,ms,'o','MarkerEdgeColor',[col3(3,:)]);
    
box on

%  ylabel('PSD $\tilde{S}^{\rm ol}_{\rm yy}$ (pm$^2$ Hz$^{-1}$)','Interpreter', 'Latex')
 ylabel('PSD $\tilde{S}^{\rm ol}_{\rm yy}$ ','Interpreter', 'Latex')
xlabel('$f$ (kHz)','Interpreter', 'Latex')
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',18,'TickLength',[0.02, 0.01], 'XMinorTick','on','xlim',...
    [141,151],'ylim', [0.08,3.5],...
    'Xtick',[142 146 150 ])
 hold off

% figure() 
% 
% for j=1:size(DATAT,2)
%     datai=DATAT.(j);
% %     plot(datai.freq, datai.psd, 'DisplayName', datai.name, 'LineWidth', 2);
%     
%      p(j)=scatter(datai.freq, datai.psd,35,'o','filled');
%       hold on
%     Ft=lorentzfit(datai.freq,datai.psd);
%      if j <size(DATAT,2)
%  f(j)=plot(datai.freq, Ft, 'DisplayName', datai.name, 'LineWidth', 2, 'color', 'k');
%  end  
%     
%    
% end
% legend
%  ylabel('PSD (pm^2/Hz)')
% xlabel('f (kHz)')
% set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',18,'TickLength',[0.02, 0.01], 'XMinorTick','on','xlim', [141,151],'ylim', [0.08,4],...
%     'Xtick',[142 146 150 ])
%  hold off
%  
%  
%  figure()
%  
%  DATA=load('data_feeding_cooling_vacuum.mat');
% DATAT=struct2table(DATA);
% 
% figure('Position',[10 20 Xpix Ypix]);
% 
% axes('Position',[bx1 0 xw2 0]/Xpix + [0 by1 0 ywi]/Ypix);  
% C=parula(4);
% for j=1:size(DATAT,2)
%     datai=DATAT.(j);
%     p(j)=scatter(datai.freq, datai.psd,55,'o','filled','MarkerFaceColor',[col3(j,:)]);
%     %,'MarkerFaceColor',[C(j,:)]
%     Ft=lorentzfit(datai.freq,datai.psd);
%  hold on
%  if j <size(DATAT,2)
%  f(j)=plot(datai.freq, Ft, 'DisplayName', datai.name, 'LineWidth', 2, 'color', 'k');
%  end  
%     
% end
% 
% legend
%%
% axes('Position',[bx1+xw2+xwi+2*bxm 0 0.6*xwi 0]/Xpix + [0 by1 0 ywi]/Ypix); 
% 
% ylabel('n$_\infty$ ','Interpreter', 'Latex')
% xlabel('p$_{\rm gas}$ (mbar)','Interpreter', 'Latex')
% 
% box on
% set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','on',...
%     'xlim',[1e-8, 1e-2],'ylim', [1e1, 1e6],'Xtick',[1e-8 1e-5 1e-2],'Ytick',[1e1 1e2 1e4 1e6],...
%      'yscale', 'log','xscale', 'log')
 
% axes('Position',[bx1+xw2+xwi+2*bxm 0 0.6*xwi 0]/Xpix + [0 by1 0 ywi]/Ypix); 
% 
% [img, ~, ImageAlpha]  = imread('Fig47.png');
% %  
% %  , 'AlphaData', ImageAlpha
% 
% image(img,'XData',[200 300],'YData',[300 400], 'AlphaData', ImageAlpha);
% set(gca, 'Xtick',[ ],'XMinorTick','off','YMinorTick','off','Ytick',[ ])

 
%  
%  hold on



%%
axes('Position',[0 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix); 




xt = [bx1,bx1+xw2+bxm];
yt = [by1+ywi+bym,by1+ywi+bym];
    
 str = {'\bf a','\bf b'};


text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])

% 
% saveas(gcf,'Fig47.eps','epsc')
% saveas(gcf,'Fig47.fig')