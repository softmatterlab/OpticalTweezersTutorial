clear all
clc
close all
set(0, 'defaultFigureRenderer', 'painters')
DATA=load('data_feeding_cooling_vacuum.mat');
DATAT=struct2table(DATA);


col3=[0.00 0.45 0.74;, 216/255 82/255 24/255; 237/255 177/255 32/255; 128/255 128/255 128/255 ];



xw2 = 400;
bx1 = 130;     % extra space at the left
bx2 = 100;     % extra space at the right
bxm = 20;     % extra space at the midle
xwi = 2*xw2+0.7*xw2+bx1+bxm;    % width of the plot square


Xpix=1400;   % total width 
ywi = 400;    % length frame with function

by1 = 100;     % extra space below
by2 = 30;     % extra space up
bym = 170      %extra space at the midle
Ypix = by1+2*ywi+2*by2+bym;  % width in pixel


figure('Position',[100 -500 Xpix Ypix]);

axes('Position',[bx1 0 xw2 0]/Xpix + [0 by1+ywi+bym 0 ywi]/Ypix);  

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

 
%  LL= legend ([p(1)  p(3) p(2) p(4) f(1)],'\gamma_{\rm fb} = 0.3 kHz','\gamma_{\rm fb} = 0.9 kHz',...
%      '\gamma_{\rm fb} = 3.0 kHz','shotnoise','Lorentzian fit','Interpreter','latex','Box','off','Position',[155 280 1 1], 'boxoff')


 LL= legend ([p(1)  p(3) p(2) p(4) f(1)],'\gamma_{\rm fb} =0.3 kHz','\gamma_{\rm fb} =0.9 kHz',...
     '\gamma_{\rm fb} =3.0 kHz','shotnoise','fit','Interpreter','latex','Box','off','Position',[230 560 1 1])
 
LL.FontSize = 18;
LL.Box = 'off';
LL.NumColumns = 2;
  
 
 set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01],...
     'XMinorTick','on','xlim', [141,151],'ylim', [2e-1,4],'yscale', 'log','Ytick',[1 2 3  4])
 



 %%
 
DATA=load('data_temperature.mat')
DATAT=struct2table(DATA);
axes('Position',[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix); 

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


%  LL= legend ([p(1) p(3) p(2) ],'p = 1.4 x 10^{-8} mbar','p = 1.2 x 10^{-7} mbar',...
%      'model','Interpreter','latex','Box','off','Position',[620 120 1 1], 'boxoff')

LL= legend ([p(1) p(3) p(2) ],'p=1.4x10^{-8} mbar','p=1.2x10^{-7} mbar',...
     'model','Interpreter','latex','Box','off','Position',[825 300 1 1], 'boxoff')
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
axes('Position',[bx1+xw2+bxm 0 xw2 0]/Xpix + [0  by1+ywi+bym 0 ywi]/Ypix); 


 datai=DATAT.(1);
     Ft=lorentzfit(datai.freq,datai.psd);
   ms=55;
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
%  ylabel('PSD $\tilde{S}^{\rm ol}_{\rm yy}$ ','Interpreter', 'Latex')
xlabel('$f$ (kHz)','Interpreter', 'Latex')
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','on','xlim',...
    [141,151],'ylim', [2e-1,4],'Ytick',[],'yscale', 'log')
% 'Xtick',[142 146 150 ]
 hold off

%%
load('dataconverted.mat') ;

x={xx; yx; zx};
y={xy; yy; zy};
fd=[17 17 16];
xf=logspace(-9,-2,100)';


axes('Position',[2*bx1+2*xw2+bxm 0 0.7*xw2 0]/Xpix + [0  by1+ywi+bym 0 ywi]/Ypix); 

for i=1: 3
   
    modelFun =  @(p,x) p(1)+p(2).*x;

    startingVals = [80 1];
    coefEsts = nlinfit(x{i,1}(1:fd(i)), y{i,1}(1:fd(i)), modelFun, startingVals);


    a=coefEsts(1);
    b=coefEsts(2);
   
    
 
     p(4)=plot(xf,a+b*xf, 'LineWidth', 2, 'color', 'k')
     hold on  
     p(i)=scatter(x{i,1},y{i,1},'o','MarkerEdgeColor',[col3(i,:)])
    
    
end


 LL= legend ([p(1) p(2) p(3) p(4) ],'x','y', 'z','fit','Interpreter','latex','Box','off','Position',[960 570 1 1])
 
LL.FontSize = 18;
LL.Box = 'off';



ylabel('n$_\infty$ ','Interpreter', 'Latex')


xlabel('$P_{\rm gas}$ (mbar)','Interpreter', 'Latex')
box on
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','on',...
    'xlim',[6e-9, 1e-2],'ylim', [1e1, 1e6],'Xtick',[1e-8 1e-5 1e-2],'Ytick',[1e1 1e2 1e4 1e6],...
     'yscale', 'log','xscale', 'log')
 

% %%
% clear all; close all;
% load('dataconverted.mat') ;
% 
% x=zx;
% y=zy;
% 
% xt=(x(1:16));
% yt=(y(1:16));
% 
% 
% figure
% 
% scatter (xt,yt)
% 
% 
% 
% 
% modelFun =  @(p,x) p(1)+p(2).*x;
% 
% 
% startingVals = [80 1];
% coefEsts = nlinfit(xt, yt, modelFun, startingVals);
% 
% 
% a=coefEsts(1);
% b=coefEsts(2);
% 
% figure(1)
% 
% scatter (x,y)
% % loglog (xt,yt)
% hold on
% plot (x,a+(b*x))
% 
% set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','on',...
%     'xlim',[6e-9, 1e-2],'ylim', [1e1, 1e6],'Xtick',[1e-8 1e-5 1e-2],'Ytick',[1e1 1e2 1e4 1e6],...
%      'yscale', 'log','xscale', 'log')
%  



%%
axes('Position',[0 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix); 


xt = [bx1,bx1+xw2+bxm,2*bx1+2*xw2+bxm,bx1];
yt = [by1+2*ywi+bym+by2,by1+2*ywi+bym+by2,by1+2*ywi+bym+by2,by1+ywi+by2];
    
 str = {'\bf a','\bf b','\bf c','\bf d'};


text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])


saveas(gcf,'Fig47.eps','epsc')
saveas(gcf,'Fig47.fig')