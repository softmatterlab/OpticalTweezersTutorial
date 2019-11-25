%% Program that reads the result files and make the plot of stiffness respect to force.

clear all;close all;clc;

path='F:\Document\GitHub\tweezers_AOP_tutorial\Marta\Results'; %where Stiffness result files are located


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

figure(1)
set(gcf,'name','Stiffness')
errorbar(F1,S1,ES1,'or')
hold on;
errorbar(F2,S2,ES2,'xr')
hold on;
errorbar(F3,S3,ES3,'ob')
hold on;
errorbar(F4,S4,ES4,'xb')
set(gca,'fontsize',14)
xlabel('Force (pN)','fontsize', 18);
ylabel('Stiffness (pN/nm)','fontsize', 18);



%% Program that reads the result files and make the plot of force - extension curve.


% clear all;close all;clc;

path='F:\Document\GitHub\tweezers_AOP_tutorial\Marta\Results'; %where individual trajectory files are located


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


figure(2)
set(gcf,'name','Force-Extension')
plot(F1,E1,'r')
hold on;
plot(F2,E2,'m')
hold on;
plot(F3,E3,'b')
hold on;
plot(F4,E4,'c')
set(gca,'fontsize',14)
ylabel('Force (pN)','fontsize', 18);
xlabel('Extension (nm)','fontsize', 18);
%% 

bx1 = 100;     % bordo a sinistra
xwi = 560;    % larghezza riquadro con funzione
bx2 = 30;     % bordino a destra

Xpix = 3*bx1+3*xwi+2*bx2;  % larghezza figura in pixel
Xpix =1400;
by1 = 110;     % bordo in basso
ywi = 500;    % altezza riquadro con funzione
by2 = 50;     % bordo in alto

Ypix = by1+ywi+by2;  % larghezza figura in pixel


col3=[0.00,0.45,0.74];
col4=[0.8500, 0.3250, 0.0980];


figure('Position',[10 20 Xpix Ypix]); % crea la figura
axes('Position',[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  % fa in modo di centrare il riquadro degli assi nella posizione voluta


% errorbar(F1,S1,ES1,'o','markeredgecolor','k', 'MarkerSize',12)
errorbar(F1,S1,ES1,'.','MarkerSize',30, 'Color','r')
hold on;
errorbar(F2,S2,ES2,'^','MarkerSize',10, 'Color','r','MarkerFaceColor','r')
hold on;
errorbar(F3,S3,ES3,'.','MarkerSize',30, 'Color',col3)
hold on;
errorbar(F4,S4,ES4,'^','MarkerSize',10, 'Color',col3,'MarkerFaceColor',col3)

xlabel('Force (pN)','Interpreter','Latex', 'FontSize',30);
ylabel('Stiffness (pN/nm)','Interpreter','Latex', 'FontSize',30);

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','on');
