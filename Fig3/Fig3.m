close all, clear all;

%%  ==============Parameter declaration============


kB=1.38e-23; % Boltzmann constant [m^2kg/s^2K]
T=300;  % Temperature [K]
r=1.03E-6;      % Particle radius [m]
v=0.00002414*10^(247.8/(-140+T));  % Water viscosity [Pa*s]
gamma=pi*6*r*v; %[m*Pa*s]


%%  ==============Selecting file============
% 
% [filename,pathname,d] = uigetfile('*.mat;');
% [pathstr,name,ext] = fileparts(filename);


%%  =========Loading selected file============

X=load('X.mat');
Y=load('Y.mat');
Z=load('Z.mat');


%%  ==============Reading data===============
N=length(X.Vx);

nd=3;
nbins=60;
% dt=A(2,1)-A(1,1);
dt=1E-5;
fs=1/dt;
t=[0:dt:dt*(N-1)]';

    
% t=A(:,1);
%     x=A(:,2); 
    x=X.Vx(:,nd); 
    x=x-mean(x);

   
%     y=A(:,3);
    y=Y.Vy(:,nd);
    y=y-mean(y);
    
%     z=A(:,4); 
     z=Z.Vz(:,nd); 
    z=z-mean(z);
      
      
      




%%  ==============X Potential analysis===============
 

% Cal=Acf;
% xmic=x/Cal.bx;
[a,b]=hist(x, nbins);
f = fit(b.',a.','gauss1');
cfx = coeffvalues(f);
a1=cfx(1);
b1=cfx(2);
c1=cfx(3);

fnx=  a1*exp(-((b-b1)/c1).^2);


    ptx=-log(a);

         
    modelFun =  @(p,x) 0.5*((x+p(1)).^2).*p(2)+p(3);

startingVals = [0 10 11.07];
coefEsts = nlinfit(b, ptx, modelFun, startingVals);

Pot.x_0=coefEsts(1);
Pot.k_x=coefEsts(2);
Pot.Ax_0=coefEsts(3);
bl=[-0.27:0.005:0.265];
cppx=0.5*((bl+Pot.x_0).^2).*Pot.k_x+Pot.Ax_0;


Pot.kx=Pot.k_x*1e12*kB*T/1E6;



%%

col3=[0.00,0.45,0.74];

bx1 = 130;     % extra space at the left
xwi = 335;    % width of the plot square
bx2 = 100;     % extra space at the right

Xpix=1400;   % total width in pixel

by1 = 100;     % extra space below
ywi = 300;    % high frame with function
by2 = 30;     % extra space up

Ypix = by1+ywi+2*by2;  % total high in pixel




figure('Position',[2000 20 Xpix Ypix]); % crea la figura
axes('Position',[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  
hold on;
plot(t,x*1000,'color', col3,'LineWidth',  1);  



box on

xlabel('$t\ (\rm s)$','Interpreter','Latex','FontSize',30 );
ylabel('$x\ (\rm{nm})$','Interpreter','Latex','FontSize',30);

%  xticks(0:1:10);
%    ylim([-300 300]);
%    yticks(-300:100:300);
 set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'Xtick',[0:2:10],'Ytick',[-300:100:300],'ylim', [-300 300]);
%  xlim([-0.3 0.3001]);




%%
axes('Position',[(2*bx1+xwi) 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
hold on;



%  scatter(b(5:(end-2))*1000,ptx(5:(end-2)),80,'o')


 scatter(b*1000,ptx,80,'o','markerfacecolor', col3,'LineWidth',0.5)

% set(gca,'xtick', 0.1)


hold on
% plot(b(5:(end-2)),cppx(5:(end-2)),'LineWidth',3,'Color',col2)
plot(bl*1000,cppx,'LineWidth',3,'Color','r')

% 
% xt = [-0.4];
% yt = [3 ];
% str = {'c'};
% text(xt,yt,str,'FontSize',34)



box on
%   xlim([-0.3 0.3001])
%      set(gca,'XTick',[Min : 0.1 : Max]);
 xticks(-300:150:300);
  xlim([-300 300]);
  yticks(-12:2:2);
   ylim([-12 2 ]);


xlabel('$x\ (\rm nm)$','Interpreter','Latex')

ylabel('$U\ (k_{\rm B}T)$','Interpreter','Latex')

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5, 'FontSize',25);

 %% 
axes('Position',[(3*bx1+2*xwi-30) 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
hold on;
      

 bar(b*1000,a)
 hold on 
 plot(b*1000,fnx,'LineWidth',3,'Color','r')
 box on
 xticks(-300:150:300);
 xlim([-300 300]);
 yticks(-0:1E4: 7.0001E4);
 ylim([-0 7E4]);


ylabel('$\rho\ (\rm {counts})$','Interpreter','Latex')
xlabel('$x\ (\rm nm)$','Interpreter','Latex')

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25);




%%
 axes('Position',[0 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
hold on

xt = [bx1,2*bx1+xwi,3*bx1+2*xwi-95];
yt = [ by1+ywi+by2,by1+ywi+by2,by1+ywi+by2];
str = {'\bf a','\bf b','\bf c'};
text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])

%%

saveas(gcf,'Fig3.eps','epsc')
saveas(gcf,'Fig3.fig')