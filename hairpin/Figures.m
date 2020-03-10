%% Program that reads the result files and make the plot of stiffness respect to force.

clear all;close all;clc;

path='Results'; %where Stiffness result files are located


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


clear all;close all;clc;

path='Results'; %where individual trajectory files are located


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

