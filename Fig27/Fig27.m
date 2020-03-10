clear all
close all
load('data_doblePozo_CP.mat');

ds=0.02e-6;% Defines the grid in m
dd=0.07e-6;% defines the radius of the domain to estimate the forces.
Nlim=5000;

gxv=min_x+ds/2:ds:max_x-ds/2;
gyv=min_y+ds/2:ds:max_y-ds/2;
addpath 'general'
addpath '/home/alejandro/Documents/MATLAB/general/'
nfile='figsDoubleWell';
mkdir(nfile);
%Cajas
[GX,GY]=meshgrid(gxv,gyv);
lx=length(gxv);
ly=length(gyv);
boxes=[]; 
nb=1;
for j=1:ly
    disp([nb,(lx)*(ly)])
    for k=1:lx
        % buscamos todos los datos que pertenecen a la cajita j
        XY0=[GX(j,k),GY(j,k)];
        [kx_mle(j,k),ky_mle(j,k),omega1(j,k),omega2(j,k),Difx(j,k),Dify(j,k),Nd(j,k),Fx(j,k),Fy(j,k)]=...
            mlemethod_v3_full(Vx',Vy',dt,T,XY0,dd,Nlim);
        nb=nb+1;
    end
end
save('MapafuerzasDoblePozo.mat')


%send Dx, Dy 
%%
% 


addpath 'general'


load('data_doblePozo_CP.mat');
KB=1.38064852e-23;%m^kgs^-2K-1

ds=0.02e-6;% Defines the grid in m
dd=0.07e-6;%0.12e-6;% defines the domain to estimate the forces.
Nlim=5000;
    
gxv=min_x+ds/2:ds:max_x-ds/2;
gyv=min_y+ds/2:ds:max_y-ds/2;


r1=[1.589e-6,1.959e-6];%A
r2=[1.979e-6,1.929e-6];%S
r3=[2.294e-6, 1.944e-6]; %B
rc=[r1;r2;r3];


y0=1.9433e-6;

%l
lx=length(gxv);


Xs=Vx;
Ys=Vy;


% Gets the equilibrum points
for k=1:3
    XY0=[rc(k,1),y0];
    [kx_mle(k),ky_mle(k),omega1(k),omega2(k),Difx(k),Dify(k),Nd(k),Fx(k),Fy(k)]=...
        mlemethod_v3_full(Xs',Ys',dt,T,XY0,dd,Nlim);
end


%starts the potential method

Xedges = gxv;
Yedges = gyv;
figure(4)
hh=histogram2(Xs,Ys,Xedges,Yedges,'Normalization','probability');
%potxy=-Kb*T*log(hh.Values(ceil(length(Yedges)/2),:));
Nh=hh.Values;
in=find(Nh==0);
Nh(in)=NaN;
potxy=-KB*T*log(Nh);
xp=(gxv(2:end)+gxv(1:end-1))/2;
yp=(gyv(2:end)+gyv(1:end-1))/2;
[X,Y]=meshgrid(xp,yp);

DUA=(interp2(X,Y,potxy',rc(2,1),rc(2,2),'linear')-interp2(X,Y,potxy',rc(3,1),rc(3,2),'linear'))/(KB*T);

disp('DUA')
disp(DUA);
DUB=(interp2(X,Y,potxy',rc(2,1),rc(2,2),'linear')-interp2(X,Y,potxy',rc(1,1),rc(1,2),'linear'))/(KB*T);
disp('DUB')
disp(DUB);
potxy=(potxy-interp2(X,Y,potxy',rc(2,1),rc(2,2),'linear'))/(KB*T);



%[X0,Y0]=meshgrid(y0,yp);
[X0,Y0]=meshgrid(xp, y0);
% [X10,Y10]=meshgrid(xp,rc(1,2));
% [X20,Y20]=meshgrid(xp,rc(2,2));
% [X30,Y30]=meshgrid(xp,rc(3,2));

[X10,Y10]=meshgrid(rc(1,1), yp);
[X20,Y20]=meshgrid(rc(2,1), yp);
[X30,Y30]=meshgrid(rc(3,1),yp);

poty0=interp2(X,Y,potxy',X0,Y0,'linear');

potx1=interp2(X,Y,potxy',X10,Y10,'linear');
potx2=interp2(X,Y,potxy',X20,Y20,'linear');
potx3=interp2(X,Y,potxy',X30,Y30,'linear');


save('PerfilPotencialFuerzaDoblePozo.mat','X0','Y0','poty0',...
    'X10','Y10','potx1',...
    'X20','Y20','potx2',...
    'X30','Y30','potx3',...
    'kx_mle','ky_mle','Fx','Fy','y0','rc','DUA','DUB')

%%
close all
clear all
xwis = 230;    % width of the plot square
bx1 = 100;     % extra space at the left
bx2 = 15;     % extra space at the right
xwil=450;
xwit=3*xwis+2*bx2;

Xpix = xwit+2*bx1+2*bx2+xwil;  % total

ywi = 210;    % length of the plot  
by1 = 70;     % extra space below
by2 = 25;     % extra space up

Ypix = 3*by1+3*ywi+3*by2;



% figure of the quiver and probability distribution
figure('Position',[10 20 Xpix Ypix]);
load('PerfilPotencialFuerzaDoblePozo.mat')
load('MapafuerzasDoblePozo.mat')
addpath 'general'
%quiver and prob

axes('Position',[bx1 0 xwit 0]/Xpix + [0 3*by1+2*ywi+2*by2 0 ywi]/Ypix);
text(-10,-10,'a','Interpreter','Latex','FontSize',34)
%set(gcf,'Position',[600 300 600 900])
%axes('Position',[0.16 0.7 1-0.17 0.5-0.18])
imagesc(gxv*1e6,gyv*1e6,Nd);

hold on
quiver(GX'*1e6,GY'*1e6,Fx',Fy','r')
plot(rc(3,1)*1e6,rc(3,2)*1e6,'ok','MarkerFaceColor','k')
plot(rc(2,1)*1e6,rc(2,2)*1e6,'ok','MarkerFaceColor','White')
plot(rc(1,1)*1e6,rc(1,2)*1e6,'ok','MarkerFaceColor','k')

ylabel('$x(\mu \rm m)$','Interpreter','Latex','FontSize',30)
xlabel('$y(\mu \rm m)$','Interpreter','Latex','FontSize',30)
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
axis([1.3 2.5 1.7 2.2])
axis tight
hold off

%%
load('PerfilPotencialFuerzaDoblePozo.mat')

KB=1.38064852e-23;%m^kgs^-2K-1
T=273.15+22;
%x-potential
axes('Position',[bx1 0 xwit 0]/Xpix + [0 2*by1+ywi+by2 0 ywi]/Ypix);
plot(X0*1e6,poty0,'LineWidth',2)
potf=@(x)spline(X0,poty0,x);
hold on
% plot(rc(1,2)*1e6,potf(rc(1,2)),'ok','MarkerFaceColor','k')
% plot(rc(2,2)*1e6,potf(rc(2,2)),'ok','MarkerFaceColor','White')
% plot(rc(3,2)*1e6,potf(rc(3,2)),'ok','MarkerFaceColor','k')
plot(rc(1,1)*1e6,potf(rc(1,1)),'ok','MarkerFaceColor','k')
plot(rc(2,1)*1e6,potf(rc(2,1)),'ok','MarkerFaceColor','White')
plot(rc(3,1)*1e6,potf(rc(3,1)),'ok','MarkerFaceColor','k')
xx=(-0.2:0.01:0.2)*1e-6;
plot((xx+rc(1,1))*1e6,1/2*kx_mle(1)*xx.^2/(KB*T)+potf(rc(1,1)),'--k','LineWidth',2)
plot((xx+rc(2,1))*1e6,1/2*kx_mle(2)*xx.^2/(KB*T)+potf(rc(2,1)),'--k','LineWidth',2)
plot((xx+rc(3,1))*1e6,1/2*kx_mle(3)*xx.^2/(KB*T)+potf(rc(3,1)),'--k','LineWidth',2)

text(rc(3,1)*1e6-0.02,potf(rc(3,1))+1.5,'A','Interpreter','Latex','FontSize',34, 'Color' ,'r')
text(rc(2,1)*1e6-0.02,potf(rc(2,1))+1.5,'S','Interpreter','Latex','FontSize',34, 'Color' ,'r')
text(rc(1,1)*1e6-0.02,potf(rc(1,1))+1.5,'B','Interpreter','Latex','FontSize',34, 'Color' ,'r')

xlabel('$$x(\mu \rm m)$$','Interpreter','Latex','FontSize',30)
ylabel('$U(x)(k_{\rm B} T)$','Interpreter','Latex','FontSize',30)
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
axis([1.3 2.5 -4 6])

%y-potential 1
axes('Position',[bx1 0 xwis 0]/Xpix + [0 by1 0 ywi]/Ypix);

xc=rc(3,2);
kk=ky_mle(3);
potx=potx3;
Y0=Y30;

plot(Y0*1e6,potx,'LineWidth',2)
hold on

potfx=@(x)spline(Y0,potx,x);
xx=(-0.12:0.01:0.12)*1e-6;
plot((xx+xc)*1e6-0.005,1/2*kk*xx.^2/(KB*T)+potfx(xc),'--k','LineWidth',2)
plot(xc*1e6-0.005,potfx(xc),'ok','MarkerFaceColor','k')
ylabel('$U(y)(k_{\rm B} T)$','Interpreter','Latex','FontSize',30)
xlabel('$y( \mu \rm m)$','Interpreter','Latex','FontSize',30)
%axis([1.75,2.15,-4,7])
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
hold off


%y-potential 2

axes('Position',[bx1+xwis+bx2 0 xwis 0]/Xpix + [0 by1 0 ywi]/Ypix);
xc=rc(2,2);
kk=ky_mle(2);
potx=potx2;
Y0=Y20;

plot(Y0*1e6,potx,'LineWidth',2)
hold on

potfx=@(x)spline(Y0,potx,x);
xx=(-0.12:0.01:0.12)*1e-6;
plot((xx+xc)*1e6+0.005,1/2*kk*xx.^2/(KB*T)+potfx(xc),'--k','LineWidth',2)
plot(xc*1e6+0.005,potfx(xc),'ok','MarkerFaceColor','White')

xlabel('$y(\mu \rm m)$','Interpreter','Latex','FontSize',30)
%axis([1.75,2.15,-4,7])
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
yticks([])
hold off
%y-potential 3
axes('Position',[bx1+2*xwis+2*bx2 0 xwis 0]/Xpix + [0 by1 0 ywi]/Ypix);
xc=rc(1,2);
kk=ky_mle(1);
potx=potx1;
Y0=Y10;

plot(Y0*1e6,potx, 'LineWidth',2)
xlabel('$y(\mu \rm m)$','Interpreter','Latex','FontSize',30)
hold on

potfx=@(x)spline(Y0,potx,x);
xx=(-0.12:0.01:0.12)*1e-6;
plot((xx+xc)*1e6-0.005,1/2*kk*xx.^2/(KB*T)+potfx(xc),'--k','LineWidth',2)
plot(xc*1e6-0.005,potfx(xc),'ok','MarkerFaceColor','k')
axis([1.75,2.15,-4,7])
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
yticks([])

hold off


%%



load('data_doblePozo_CP.mat');

addpath general
%(genpath('/home/alejandro/Documents/MATLAB/calPot'))

% Ubicacion de puntos criticos.



r1=[1.589e-6,1.959e-6];%A
r2=[1.979e-6,1.929e-6];%S
r3=[2.294e-6, 1.944e-6]; %B
rc=[r1;r2;r3];
%ds=0.02e-6;% Defines the grid in m
dd=0.07e-6;% defines the radius of the domain to estimate the forces.

% FInd the stiffness in the critical points

for j=1:3
    XY0=[rc(j,1),rc(j,2)];
    [kx_mle(j),ky_mle(j),omega1(j),omega2(j),Difx(j),Dify(j),Nd(j),Fx(j),Fy(j)]=...
            mlemethod_v3_full(Vx',Vy',dt,T,XY0,dd,5000);
end

gamma=6*pi*eta*a;
tau0_A=pi*gamma*sqrt(ky_mle(2))/(sqrt(abs(kx_mle(2))*ky_mle(1)*kx_mle(1)));
tau0_B=pi*gamma*sqrt(ky_mle(2))/(sqrt(abs(kx_mle(2))*ky_mle(3)*kx_mle(3)));
disp(' inverse of relaxation times of each well')
disp('tau_0')
%disp([tau0_A,tau0_B])
disp([1/tau0_A,1/tau0_B])   
tau01d_A=pi*gamma/sqrt(abs(kx_mle(2))*kx_mle(1));
tau01d_B=pi*gamma/sqrt(abs(kx_mle(2))*kx_mle(3));
%disp([tau01d_A,tau01d_B])


save('CriticalPoints.mat')
%% Compute the mean passage or residence time



load('CriticalPoints.mat','tau0_A','tau0_B','tau01d_A','tau01d_B')
load('PerfilPotencialFuerzaDoblePozo.mat','DUA','DUB')

mtauA=tau0_A*exp(DUA);
mtauB=tau0_B*exp(DUB);
disp('inverse theoretical mean residence time')
disp('1/tauA_m, 1/tauB_m')
disp([1/mtauA,1/mtauB])

%mtau1dA=tau01d_A*exp(DUA);
%mtau1dB=tau01d_B*exp(DUB);

%disp([mtau1dA,mtau1dB])


%% Computes the probability distribution



load('CriticalPoints.mat')
addpath general
%time series
axes('Position',[bx1+xwit+bx2+bx1 0 xwil 0]/Xpix + [0 3*by1+2*ywi+2*by2 0 ywi]/Ypix);



x0=Vx(1:end-1)-rc(2,1);
y0=Vy(1:end-1)-rc(2,2);
 
vect=(0:(length(x0)-1))*dt;

plot(vect,(x0+rc(2,1))*1e6,'k','LineWidth',1)
hold on
plot(vect,ones(size(vect))*rc(1,1)*1e6,'--k')
plot(vect,ones(size(vect))*rc(2,1)*1e6,'--k')
plot(vect,ones(size(vect))*rc(3,1)*1e6,'--k')
text(0,rc(3,1)*1e6+0.04,'A','Interpreter','Latex','FontSize',34, 'Color' ,'r')
text(0,rc(2,1)*1e6+0.04,'S','Interpreter','Latex','FontSize',34, 'Color' ,'r')
text(0,rc(1,1)*1e6+0.04,'B','Interpreter','Latex','FontSize',34, 'Color' ,'r')

axis([0,40,1.3,2.6])
xlabel('$t(s)$','Interpreter','Latex','FontSize',30)
ylabel('$x(\mu \rm m)$','Interpreter','Latex','FontSize',30)
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
%set(gcf,'Position',[676 341 909 622])

% simple cross
% Find the equilibrium points
indpos=find(x0>0);
x0_pos=mean(x0(indpos));
s0_pos=std(x0(indpos));
indneg=find(x0<0);
x0_neg=mean(x0(indneg));
s0_neg=std(x0(indneg));

% Conditions for crossing: cross zero, cross from an equilibrium point,
% this mean that the points before crossing have to have a minimum average
% near w_(pos,neg). The crossing is defined by the mean of the average. 

n_pos=0;
n_neg=0;
kcross=find(x0(1:end-1).*x0(2:end)<0);
for j=1:length(kcross)-1
    my=mean(x0(kcross(j):kcross(j+1)));    
    if my>x0_pos-s0_pos
        % yes, it is a cross event from the well in the positive side
        n_pos=n_pos+1;
        tau_pos(n_pos)=vect(kcross(j+1))-vect(kcross(j));
        k_pos(n_pos)=kcross(j+1);
    else
        if my<x0_neg+s0_neg
            % yes, it is a crossing event from the negative side
          n_neg=n_neg+1;
         tau_neg(n_neg)=vect(kcross(j+1))-vect(kcross(j));
          k_neg(n_neg)=kcross(j+1);

        end
    end

end

hold off

% distance between wells
dw=rc(3,1)-rc(1,1);
disp('distance')
disp(dw)
%%
%close all
%dist A
axes('Position',[bx1+xwit+bx2+bx1 0 xwil 0]/Xpix + [0 2*by1+ywi+by2 0 ywi]/Ypix);
h2=histogram(tau_neg,20,'Normalization','pdf');
tt=linspace(0,h2.BinEdges(end),100);
mtau_neg=mean(tau_neg);
hold on
plot(tt,1/mtau_neg*exp(-tt/mtau_neg),'LineWidth',2);

xlabel('$\tau_{A}(s)$','Interpreter','Latex','FontSize',30)
ylabel('$\rho_A$','Interpreter','Latex','FontSize',30)
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);

axis([0,22,0,0.35])
%dist B
axes('Position',[bx1+xwit+bx2+bx1 0 xwil 0]/Xpix + [0 by1 0 ywi]/Ypix);
h1=histogram(tau_pos,20,'Normalization','pdf');
tt=linspace(0,h1.BinEdges(end),100);
mtau_pos=mean(tau_pos);
hold on
plot(tt,1/mtau_pos*exp(-tt/mtau_pos),'LineWidth',2);
xlabel('$\tau_{B}(s)$','Interpreter','Latex','FontSize',30)
ylabel('$\rho_B$','Interpreter','Latex','FontSize',30)
set(gca,'TickLabelInterpreter','tex', 'linewidth',1.5, 'FontSize',25);
axis([0,4.5,0,1.5])
%experimental mean residence time
disp('inverse experimental mean residence time')
disp([1/mtau_pos, 1/mtau_neg])

axes('Position',[(0) 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
hold on
xwi=xwit;
xt = [bx1-85,bx1+xwi+bx2,bx1-85,bx1+xwi+bx2,bx1-85,bx1+xwi+bx2];
yt = [ 3*by1+3*ywi+2*by2, 3*by1+3*ywi+2*by2,2*by1+2*ywi+2*by2 , 2*by1+2*ywi+2*by2,by1+ywi+by2, by1+ywi+by2];
str = {'\bf a','\bf b','\bf c', '\bf d','\bf e','\bf f'};
text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])