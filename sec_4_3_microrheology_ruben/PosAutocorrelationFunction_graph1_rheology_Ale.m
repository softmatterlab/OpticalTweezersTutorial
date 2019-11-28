%COMPUTE velocity autocorrelation  of particle trapped by optical tweezers
%simple fluid and determine viscosity from decay time

% Files Water.txt, PNP.txt and CPyCl4mM.txt are needed to run PosAutocorrelationFunction.m and get the plots and fits 
% of Figure 2(a) and (b) on the manuscript. The first column of these files contains the particle position (a), while the 
% second has the correspondinf times. A .mat version of these data is also given with the same information.

clear all
close all

T0=22; %bath temperature
kB=1.38e-23; %Boltzman constant
T0 =T0 + 273.16;
R = [2.73, 3.25, 3.25]*1e-6/2; %particle radius
R = 2.73*1e-6/2;
% 
xwi = 540;    % width of the plot square
bx1 = 140;     % extra space at the left
bx2 = 20;     % extra space at the right

Xpix = 2*xwi+2*bx1+2*bx2;  % total
Xpix=1400;
ywi = 150;    % length riquadro con funzione
YYWi=3*ywi;
by1 = 95;     % extra space below
by2 = 50;     % extra space up

Ypix = 1*by1+1*YYWi+1*by2;  % larghezza figura in pixel
%number of bins of the histogram, if not set default is 50
figure('Position',[10 20 Xpix Ypix]);



  SS=[{'o'},{'s'},{'d'}];


Extension='.txt';
Filepath = 'PassiveMicrorheologyData/';
Filename{1}=['Water'];
Filename{2}= [ 'PNP'];
Filename{3}=['CPyCl4mM'];

color{1}=[1,0,0];
color{2}=[0,0.8,1];
color{3}=[0.9,0.4,0.1];
for jj=1:3
filname = [Filepath Filename{jj} Extension];
Data = dlmread(filname,'',1,0);
positionintefig1=[bx1 0 xwi 0]/Xpix + [0 by1+(3-jj)*ywi 0 ywi]/Ypix;
axes( 'Position',positionintefig1);
%Read data
[aaa,bbb]=size(Data);
x = Data(:,1); %position in meters
t = Data(:,2); %time in seconds
plot(t,x/1e-9, 'Color', color{jj}, 'LineWidth', 1.5)
xlim([0 10])
xticks([0, 2, 4, 6, 8 ,10]);
yticks([-200, 0 , 200]);


ylim([-300 300])


 
ylabel('$x ( \rm nm$)','Interpreter','Latex', 'FontSize',30);

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01]);
%   legend(Filename(jj))
%   legend boxoff 
 pos=[0.62 0.36 0.1];
  LL= legend ({convertStringsToChars(Filename{jj})},'Interpreter','latex','Box','off','Fontsize', 18,'Position',[0.08 pos(jj) 0.5 0.5])
  
% 
if (jj==1 || jj==2)
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5, 'FontSize',25,'Xticklabel',[]);

end 
end

%  
%  LL= legend ('$\rm PNP$','Interpreter','latex','Box','off','Position',[0.87 0.38 0.1 0.2])
%  LL.FontSize = 18
%  
%  LL= legend ('$\rm CPyCl4mM$','Interpreter','latex','Box','off','Position',[0.97 0.28 0.1 0.2])
%  LL.FontSize = 18
xlabel('$t\rm (s)$', 'Interpreter','Latex', 'FontSize',30)
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01]);
%title('Stochastic trajectory of trapped particle')
%%
positionintefig2=[2*bx1+xwi+bx2 0 xwi 0]/Xpix + [0 by1 0 YYWi]/Ypix;
axes( 'Position',positionintefig2);

%Compute autocorrelation function
for jj=1:3
filname = [Filepath Filename{jj} Extension];
Data = dlmread(filname,'',1,0);

%Read data
[aaa,bbb]=size(Data);
x = Data(:,1); %position in meters
t = Data(:,2); %time in seconds
xx=xcorr(x,'unbiased');

aa = length(x);

kk=100; %take only the kk-nth part of total time interval of autocorrelation function

xxcorr=xx(aa:end); %take possitive values of time for autocorrelation function 

%figure(2)

plot(t(1:floor(aa/kk)),xxcorr(1:floor(aa/kk))/1e-18,SS{jj},  'MarkerSize',8, 'Color', color{jj},'MarkerFaceColor', color{jj})
hold on
%title('Stochastic trajectory of trapped particle')

%Compute trap stiffness
varx = var(x); %variance of particle position
k = kB*T0/varx %trap stiffness from equipartition 


%fit first points of autocorrelatio function to exponential and determine
%viscosity

tt = t(1:floor(floor(aa/kk)/10));
logxxcorr = log(xxcorr(1:floor(floor(aa/kk)/10)));

coeff = polyfit(tt(1:10),logxxcorr(1:10),1)

eta = k/6/pi/R/(-coeff(1))  %viscosity in Pa*s
xxcorrfit = exp(coeff(2))*exp(coeff(1)*t(1:floor(aa/kk)));

%figure(2)
plot(t(1:floor(aa/kk)),xxcorrfit/1e-18,'--', 'LineWidth',2,'Color','k','HandleVisibility','off')

xlim([0 0.75])
end

xlabel('$t\rm (s)$', 'Interpreter','Latex', 'FontSize',30)
ylabel('$\langle x(t)x(0)\rangle(\rm nm ^2)$','Interpreter','Latex', 'FontSize',30);
%set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25,'Xticklabel',[]);
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01]);
%%

positionintefig3=[2*bx1+xwi+220 0 xwi/2+40 0]/Xpix + [0 by1/2+2*YYWi/3-50 0 YYWi/3+60]/Ypix;
axes( 'Position',positionintefig3);
%Compute autocorrelation function
for jj=1:3
filname = [Filepath Filename{jj} Extension];
Data = dlmread(filname,'',1,0);

%Read data
[aaa,bbb]=size(Data);
x = Data(:,1); %position in meters
t = Data(:,2); %time in seconds
xx=xcorr(x,'unbiased');

aa = length(x);

kk=100; %take only the kk-nth part of total time interval of autocorrelation function

xxcorr=xx(aa:end); %take possitive values of time for autocorrelation function 

%figure(2)

semilogy(t(1:floor(aa/kk)),xxcorr(1:floor(aa/kk))/xxcorr(1),SS{jj},  'MarkerSize',8, 'Color', color{jj},'MarkerFaceColor', color{jj})
hold on
%title('Stochastic trajectory of trapped particle')

%Compute trap stiffness
varx = var(x); %variance of particle position
k = kB*T0/varx %trap stiffness from equipartition 


%fit first points of autocorrelatio function to exponential and determine
%viscosity

tt = t(1:floor(floor(aa/kk)/10));
logxxcorr = log(xxcorr(1:floor(floor(aa/kk)/10)));

      coeff = polyfit(tt,logxxcorr,1);

eta = k/6/pi/R/(-coeff(1));  %viscosity in Pa*s
xxcorrfit = exp(coeff(2))*exp(coeff(1)*t(1:floor(aa/kk)));

%figure(2)


 if jj == 1
      coeff = polyfit(tt(1:10),logxxcorr(1:10),1);
 else
  coeff = polyfit(tt,logxxcorr,1);
     end

eta = k/6/pi/R/(-coeff(1));  %viscosity in Pa*s
xxcorrfit = exp(coeff(2))*exp(coeff(1)*t(1:floor(aa/kk)));
semilogy(t(1:floor(aa/kk)),xxcorrfit/xxcorrfit(1),'--', 'LineWidth',2,'Color','k')

xlim([0 0.4])
ylim([3e-2 1])
end
%legend
xlabel('$t\rm (s)$', 'Interpreter','Latex', 'FontSize',30)
ylabel('$k \langle x(t)x(0)\rangle/ k_{\rm B} T$','Interpreter','Latex', 'FontSize',30);
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',18,'TickLength',[0.02, 0.01]);
%set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25,'Xticklabel',[]);



axes('Position',[(0) 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
hold on

xt = [bx1-105,bx1+xwi+bx2+20];
yt = [ by1+ywi+330,by1+ywi+330];
str = {'\bf a','\bf b'};
text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])