%% Program that reads the result files and make the plot of stiffness respect to force.

clear all;close all;clc;

path='../results'; %where Stiffness result files are located


fid1=fopen(char(strcat(path,'/RBC1_Stiffness_U.txt')),'r');
fid2=fopen(char(strcat(path,'/RBC1_Stiffness_D.txt')),'r');
fid3=fopen(char(strcat(path,'/RBC2_Stiffness_U.txt')),'r');
fid4=fopen(char(strcat(path,'/RBC2_Stiffness_D.txt')),'r');
data1 = textscan(fid1,'%f\t%f\t%f\t%f\n', 'CommentStyle', '#');
F1=data1{1};
EF1=data1{2};
S1=data1{3};
ES1=data1{4};
data2 = textscan(fid2,'%f\t%f\t%f\t%f\n', 'CommentStyle', '#');
F2=data2{1};
EF2=data2{2};
S2=data2{3};
ES2=data2{4};
data3 = textscan(fid3,'%f\t%f\t%f\t%f\n', 'CommentStyle', '#');
F3=data3{1};
EF3=data3{2};
S3=data3{3};
ES3=data3{4};
data4 = textscan(fid4,'%f\t%f\t%f\t%f\n', 'CommentStyle', '#');
F4=data4{1};
EF4=data4{2};
S4=data4{3};
ES4=data4{4};

% figure(1)
% set(gcf,'name','Stiffness')
% errorbar(F1,S1,ES1,'or')
% hold on;
% errorbar(F2,S2,ES2,'xr')
% hold on;
% errorbar(F3,S3,ES3,'ob')
% hold on;
% errorbar(F4,S4,ES4,'xb')
% set(gca,'fontsize',14)
% xlabel('$f_{\rm trap}$(pN)','fontsize', 18);
% ylabel('Stiffness (pN/nm)','fontsize', 18);


 
bx1 = 100;     % bordo a sinistra
xwi = 560;    % larghezza riquadro con funzione
bx2 = 30;     % bordino a destra

  
Xpix =1400; % larghezza figura in pixel
by1 = 110;     % bordo in basso
ywi = 500;    % altezza riquadro con funzione
by2 = 50;     % bordo in alto

Ypix = by1+ywi+by2;  % larghezza figura in pixel


col3=[0.00,0.45,0.74];
col4=[0.8500, 0.3250, 0.0980];


figure('Position',[10 20 Xpix Ypix]); % crea la figura
axes('Position',[2*bx1+xwi+50 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix); 
% errorbar(F1,S1,ES1,'o','markeredgecolor','k', 'MarkerSize',12)
errorbar(F1,S1,ES1,'.','MarkerSize',30, 'Color','r')
hold on;
errorbar(F2,S2,ES2,'^','MarkerSize',10, 'Color','m','MarkerFaceColor','m')
hold on;
errorbar(F3,S3,ES3,'.','MarkerSize',30, 'Color',col3)
hold on;
errorbar(F4,S4,ES4,'^','MarkerSize',10, 'Color','k','MarkerFaceColor','k')


xlabel('$f_{\rm trap}$(pN)','Interpreter','Latex', 'FontSize',30);
ylabel('RBC Stiffness (pN nm$^{-1}$)','Interpreter','Latex', 'FontSize',30);
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','off');


LL= legend ({'RBC1 stretching','RBC1 releasing','RBC2 stretching','RBC2 releasing'},'Interpreter','latex','Box','off','Position',[0.63 0.69 0.1 0.2])
LL.FontSize = 18


%% Program that reads the result files and make the plot of force - extension curve.


% clear all;close all;clc;

path='../results'; %where individual trajectory files are located


fid1=fopen(char(strcat(path,'/RBC1_2T_U.txt')),'r');
fid2=fopen(char(strcat(path,'/RBC1_2T_D.txt')),'r');
fid3=fopen(char(strcat(path,'/RBC2_2T_U.txt')),'r');
fid4=fopen(char(strcat(path,'/RBC2_2T_D.txt')),'r');
data1 = textscan(fid1,'%f\t%f\t%f\t%f\n', 'CommentStyle', '#');
F1=data1{1};
E1=data1{2};
data2 = textscan(fid2,'%f\t%f\t%f\t%f\n', 'CommentStyle', '#');
F2=data2{1};
E2=data2{2};
data3 = textscan(fid3,'%f\t%f\t%f\t%f\n', 'CommentStyle', '#');
F3=data3{1};
E3=data3{2};
data4 = textscan(fid4,'%f\t%f\t%f\t%f\n', 'CommentStyle', '#');
F4=data4{1};
E4=data4{2};
% 
% 
% figure(2)
% set(gcf,'name','Force-Extension')
% plot(F1,E1,'r')
% hold on;
% plot(F2,E2,'m')
% hold on;
% plot(F3,E3,'b')
% hold on;
% plot(F4,E4,'c')
% set(gca,'fontsize',14)
% ylabel('Force (pN)','fontsize', 18);
% xlabel('Extension (nm)','fontsize', 18);




%% 




axes('Position',[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  % fa in modo di centrare il riquadro degli assi nella posizione voluta



plot(F1,E1,'r','linewidth',2)
hold on;
plot(F2, E2,'m','linewidth',2)
hold on;
plot(F3,E3,'Color',col3,'linewidth',2)
hold on;
plot(F4,E4,'k','linewidth',2)

ylabel('$f_{\rm trap}$(pN)','Interpreter','Latex', 'FontSize',30);
xlabel('$\Delta x_{\rm cell}$(nm)','Interpreter','Latex', 'FontSize',30);
 ylim([0 20])

 

LL= legend ({'RBC1 stretching','RBC1 releasing','RBC2 stretching','RBC2 releasing'},'Interpreter','latex','Box','off','Position',[0.13 0.69 0.1 0.2])
LL.FontSize = 18



set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','off');



axes('Position',[(0) 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
hold on

xt = [bx1-85,bx1+xwi+60];
yt = [ by1+ywi+30,by1+ywi+30];
str = {'\bf a','\bf b'};
text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])