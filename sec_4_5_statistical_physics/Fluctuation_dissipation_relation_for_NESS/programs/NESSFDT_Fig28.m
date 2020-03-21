%This program computes from the experimental trajectories (theta) as a function 
%of time (t) of a spherical particle (radius r) moving in water on a circle
%(radius a, friction coefficient gamma) in a non-equilibrium steady state (NESS) 
%under the action of constant non-conservative force and periodic
%potential, the following quantities:
       %-the probability current (J)
       %-the probability density function (pdfang) of the polar angle (ang)
       %-the non-conservative force F (normalized by gamma*a)
       %-the profile of the periodic potential U.
%Once J and pdfang are determined, it also computes, for the periodic observable Q = sin(ang):
       %-the autocorrelation function (C) of Q
       %-the corrective correlation term (B) due to the non-vanishing J in
        %a NESS
       %-the integrated linear response function Chi, and the corresponding linear
        %response function R, to a change in the amplitude of the potential
        %energy around a NESS.
%and compares them with the quantities Chi1 and R1 that would be observed in thermal equilibrium (F=0) 

clear all
close all


set(0, 'defaultFigureRenderer', 'painters')

%read data
 Filepath = ''
Filename = '../data/trajNESS_data';
Extension = '.txt';

T0=20; %bath temperature
kB=1.38e-23; %Boltzman constant
T0 =T0 + 273.16;
r = (2.0e-6)/2; %particle radius (m)
a = 4.12e-6; %radius of circular trajectory where particle is restricted to move (m)
eta = 1e-3; %water viscosity
gamma = 6*pi*eta*r; %friction coefficient
D = kB*T0/gamma ; %diffusion coefficient
Dth = D/a^2; %(polar) angular diffusion coefficient along the circle

filname = [Filepath Filename Extension]; 

%Read data
Data = dlmread(filname,'',1,0);

t = Data(:,end); %time (seconds)
theta = Data(:,1:end-1); %angular position (radians)

fs = 1/(t(2)-t(1)); %sampling frequency


%compute probaility current J
thetamean = mean(theta,2); %ensemble average of trajectories
coeff = polyfit(t,thetamean,1);
thetameanfit = polyval(coeff,t,1); %mean trajectory with drift due to probability current
J = coeff(1)/2/pi %probability current 


figure(1)
figure1 = figure(1);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
xlim(axes1,[0 66.66]);
ylim(axes1,[0 10*pi]);
plot(t,theta)
hold on
xlabel({'{\itt} [s]'});
ylabel({'\theta [rad]'});
box(axes1,'on');
set(axes1,'YTick',[0 2*pi 4*pi 6*pi 8*pi 10*pi],'YTickLabel',{'0','2\pi','4\pi','6\pi','8\pi','10\pi'});
plot(t,thetameanfit,'LineWidth',3,'LineStyle','--','Color',[0 0 0])


%Compute pdf of polar angle theta -> th defined from 0 to 2*pi
theta1 = mod(theta,2*pi);

figure(2)
figure2 = figure(2);
axes2 = axes('Parent',figure2);
hold(axes2,'on');
xlim(axes2,[0 66.66]);
ylim(axes2,[0 2*pi]);
plot(t,theta1(:,1))
xlabel({'{\itt} [s]'});
ylabel({'\theta [rad]'});
box(axes2,'on');
set(axes2,'YTick',[0 2*pi],'YTickLabel',{'0','2\pi'});

[L,N] = size(theta1); %N = number of different realizations of theta1
                      %L = number of points in each trajectory

th=[];
for j = 1:N
    th = [th; theta1(:,j)];
end

Ntot = length(th); %total number of NESS points to compute pdfang
M = 100; %number of bins
dang = 2*pi/M; %bin size (rad)
ang = zeros(1,M+1); % define bin vector
pdfang = zeros(1,M+1); %define pdf vector

ang(1) = 0; %first bin
ang(M+1)= 2*pi; %last bin
JJ = find(or(th < 0.5*dang,th >= 2*pi - 0.5*dang));
NJ = length(JJ);
pdfang(1) = M/2/pi/Ntot*NJ; 
pdfang(M+1) = pdfang(1); %first and last point must be equal due to periodicity 

clear JJ

for j = 2:M
    ang(j) = (j-1)*dang;
    JJ = find(and(th >=(j-1)*dang - 0.5*dang,th < (j-1)*dang + 0.5*dang));
    NJ = length(JJ);
    pdfang(j) = M/2/pi/Ntot*NJ;
    clear JJ
end

figure(3)
yyaxis left
bar(ang,pdfang,'FaceColor',[0 1 1])
hold on
xlabel({'{\theta} [rad]'});
ylabel({'\rho_{NESS}(\theta) [rad^{-1}]'});


%Compute non-conservative force F (rad / s, normalized by gamma*a) 
F = J*dang/2/pi*(sum(1./pdfang) - 1/pdfang(1))

%Compute profile H of periodic potential (rad^2 / s, normalized by gamma*a^2)
for k = 1:length(pdfang)
    H(k) = -Dth*log(pdfang(k)) + F*ang(k) - J*dang*(sum(1./pdfang(1:k))-0.5/pdfang(k)-0.5/pdfang(1)); 
end
U=gamma*a*a*H; %potential energy (in Joules)
figure(3)
yyaxis right
plot(ang,U/kB/T0,'LineWidth',3,'Color',[1 0 0])
ylabel({'U(\theta) [{\itk_BT}]'});

U0 = max(U/kB/T0) %potential barrier (in units of kB*T)


%Compute equilibrium distribution (F=0)
Bolzfact = exp(-U/kB/T0);
Z=dang*(sum(Bolzfact) - 0.5*(Bolzfact(1)+Bolzfact(end)));
pdfeq = Bolzfact/Z; %equilibrium distribution

figure(4)
figure4 = figure(4);
axes4 = axes('Parent',figure4);
hold(axes4,'on');
xlim(axes4,[0 2*pi]);
ylim(axes4,[0 4]);
area(ang,pdfang,'FaceColor',[0 1 1])
hold on
area(ang,pdfeq,'FaceColor',[0 0.498039215803146 0])
box(axes4,'on');
set(axes4,'XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','pi/2','pi','3\pi/2','pi'});
set(axes4,'YTick',[0 1 2 3 4],'YTickLabel',{'0','1','2','3','4'});
xlabel({'{\theta} [rad]'});
ylabel({'\rho_{NESS}(\theta) [rad^{-1}]'});



%Compute correlation functions C and B and 
Q = sin(theta1); %observable periodic in polar angle
dV = cos(theta1); %derivative of variable conjugate to perturbation

%Autocorrelation funcion of Q
for j = 1:N
    CQQ(:,j) = xcorr(Q(:,j),'unbiased'); %time average
end

C = mean(CQQ,2); %ensemble average
C = C(L:end); %autoccorelation function of Q for possitive time lags
tt = t(1:L)-t(1);
kk = 2;
figure(5)
figure5 = figure(5);
axes5 = axes('Parent',figure5);
hold(axes5,'on');
xlim(axes5,[0 33.33]);
ylim(axes5,[-0.1 0.3]);
plot(tt(1:floor(L/kk)),C(1:floor(L/kk)),'LineWidth',3,'Color',[0 0 1])
box(axes5,'on');
set(axes5,'XTick',[0 10 20 30],'XTickLabel',{'0','10','20','30'});
set(axes5,'YTick',[-0.1 0 0.1 0.2 0.3],'YTickLabel',{'-0.1','0','0.1','0.2','0.3'});
xlabel({'{\itt} [s]'});
ylabel({'{\itC}({\itt})'});

figure(6)
figure6 = figure(6);
axes6 = axes('Parent',figure6);
hold(axes6,'on');
xlim(axes6,[0 33.33]);
ylim(axes6,[0 0.35]);
plot(tt(1:floor(L/kk)),C(1)-C(1:floor(L/kk)),'LineStyle','-.','LineWidth',3,'Color',[0 0 1])
hold on


%Compute mean local velocity
v = smoothdata(J./pdfang,'rloess'); %apply robust quadratic regression to smooth data

for j = 1:N
    thj=theta1(:,j);
    vel = spline(ang,v,thj);
    vL(:,j) = vel;
    clear vel
end

X = Q.*vL;
for j = 1:N
    BXdV(:,j) = xcorr(X(:,j),dV(:,j),'unbiased');
end

b = mean(BXdV,2);
b = b(L:end);
tt = t(1:L)-t(1);
kk = 2;


for k = 1:L
   B(k) = (sum(b(1:k)) - 0.5*(b(1)+b(k)))/fs; 
end

B = B';
    
figure(6)
plot(tt(1:floor(L/kk)),B(1:floor(L/kk)),'Linestyle','--','LineWidth',3,'Color',[1 0 0])
plot(tt(1:floor(L/kk)),C(1)-C(1:floor(L/kk))-B(1:floor(L/kk)),'Linestyle','-','LineWidth',3,'Color',[0 0.498039215803146 0])
box(axes6,'on');
set(axes6,'XTick',[0 10 20 30],'XTickLabel',{'0','10','20','30'});
set(axes6,'YTick',[0 0.1 0.2 0.3 0.4],'YTickLabel',{'0','0.1','0.2','0.3','0.4'});
xlabel({'{\itt} [s]'});


%Compute response function
Chi1 = (C(1) - C)/kB/T0;  %uncorrected integrated response function
Chi = (C(1) - C - B)/kB/T0; %corrected integrated response function
R1 = fs*diff(Chi1); %uncorrected response function
R = fs*diff(Chi);  %corrected response function

figure(7)
figure7 = figure(7);
axes7 = axes('Parent',figure7);
hold(axes7,'on');
xlim(axes7,[0 33.33]);
ylim(axes7,[-1 8]);
plot(tt(1:floor(L/kk)),gamma*a*a*R1(1:floor(L/kk)),'LineStyle','-.','LineWidth',3,'Color',[0 0 1])
hold on
plot(tt(1:floor(L/kk)),gamma*a*a*R(1:floor(L/kk)),'Linestyle','-','LineWidth',3,'Color',[0 0.498039215803146 0])
box(axes7,'on');
set(axes7,'XTick',[0 10 20 30],'XTickLabel',{'0','10','20','30'});
set(axes7,'YTick',[0 2 4 6 8],'YTickLabel',{'0','2','4','6','8'});
xlabel({'{\itt} [s]'});
ylabel({'\gamma{\ita}^2 {\itR}({\itt}) [rad^{-2}]'});

%% 

xwi = 490;    % width of the plot square
bx1 = 120;     % extra space at the left
bx2 = 20;     % extra space at the right
bxm = 150;     % extra space at the midle


Xpix=1400;   % total width 
ywi = 200;    % length frame with function

by1 = 50;     % extra space below
by2 = 30;     % extra space up
bym = 80      %extra space at the midle
Ypix = by1+2*ywi+by2+bym;  % width in pixel


col3=[0.00,0.45,0.74];

figure('Position',[10 20 Xpix 2*Ypix]); % generate te figure

%%BOX 1
axes('Position',[bx1 0 xwi 0]/Xpix + [0 by1+ywi+bym 0 ywi]/Ypix);  
% figure();
%  axes1 = axes('Parent',figure1);
%  hold(axes1,'on');

plot(t,theta)
hold on
xlabel('$t$ (s)','Interpreter','latex','FontSize',30);
ylabel('$\theta(t)$ (rad)','Interpreter','latex','FontSize',30);
% box(axes,'on');
xlim([0 66.66]);
ylim([0 10*pi]);
% yticks ([0 2*pi 4*pi  6*pi 8*pi 10*pi])
% yticks ({'0','\pi','2\pi','3\pi','4\pi','5\pi','6\pi','7\pi','8\pi','9\pi','10\pi'})
% set(gca,'YTick',[0 2*pi 4*pi 6*pi 8*pi 10*pi],'YTickLabel',{'0','2\pi','4\pi','6\pi','8\pi','10\pi'});
plot(t,thetameanfit,'LineWidth',3,'LineStyle','--','Color',[0 0 0])

  set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'XMinorTick','off','Tickdir','in','TickLength',[0.02, 0.01],...
'YTick',[0 2*pi 4*pi 6*pi 8*pi 10*pi],'YTickLabel',{'$0$','$2\pi$','$4\pi$','$6\pi$','$8\pi$','$10\pi$'});

%%BOX 1 inset
axes('Position',[bx1+67 0 xwi/2.6 0]/Xpix + [0 by1+ywi+bym+153 0 ywi/5]/Ypix);

axes2 = axes('Parent',figure2);
hold(axes2,'on');
xlim(axes2,[0 66.66]);
ylim(axes2,[0 2*pi]);
plot(t,theta1(:,1))
xlabel('$t$ (s)','Interpreter','latex','FontSize',30);
ylabel('$\theta(t)$(rad)','Interpreter','latex','FontSize',30);
box(axes2,'on');
set(axes2,'YTick',[0 2*pi],'YTickLabel',{'0','2\pi'});

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',16,'XMinorTick','off','Tickdir','in','TickLength',[0.02, 0.01],...
    'XTick',[0 20 40 60])



%%BOX 2 
axes('Position',[bx1+xwi+bxm 0 xwi 0]/Xpix + [0 by1+ywi+bym 0 ywi]/Ypix);  

yyaxis left
bar(ang,pdfang,'FaceColor',[0 1 1])
hold on
xlabel('$\theta$ (rad)','Interpreter','latex','FontSize',30);

ylabel('$\rho_{\rm NESS}(\theta)$ (rad$^{-1})$', 'Interpreter','latex','FontSize',30);


%Compute non-conservative force F (rad / s, normalized by gamma*a) 
F = J*dang/2/pi*(sum(1./pdfang) - 1/pdfang(1))

%Compute profile H of periodic potential (rad^2 / s, normalized by gamma*a^2)
% for k = 1:length(pdfang)
%     H(k) = -Dth*log(pdfang(k)) + F*ang(k) - J*dang*(sum(1./pdfang(1:k))-0.5/pdfang(k)-0.5/pdfang(1)); 
% end
U=gamma*a*a*H; %potential energy (in Joules)
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'YMinorTick','off','Tickdir','in','TickLength',[0.02, 0.01],...
        'YTick',[0 0.25 0.50 0.75 1 ],'XTickLabel',{'$0$','$0.25$','$0.50$','$0.75$','$1$'})

hold on

yyaxis right
plot(ang,U/kB/T0,'LineWidth',3,'Color','r')
ylabel('$U(\theta) (k_{\rm B}T)$','Interpreter','latex','FontSize',30,'Color','k');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'YMinorTick','off','Tickdir','in','TickLength',[0.02, 0.01],...
    'YTick',[-100 -50 0 50 100 ],'YTickLabel',{'$-100$','$-50$','$0$','$50$','$100$'},...
    'XTick',[0 pi/2 pi 3/2*pi 2*pi ],'XTickLabel',{'$0$','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'})
% ,...
%     'YTick',[0 0.25 0.50 0.75 1 ],'XTickLabel',{'$0$','$0.25$','$0.50$','$0.75$','$1$'})




%%BOX 2 inset

axes('Position',[bx1+xwi+460 0 xwi/4.2 0]/Xpix + [0 by1+ywi+bym+110 0 ywi/2.6]/Ypix);


xlim([0 2*pi]);
ylim([0 4]);
area(ang,pdfang,'FaceColor',[0 1 1])
hold on
area(ang,pdfeq,'FaceColor',[0 0.498039215803146 0])
box(axes4,'on');
set(axes4,'XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','pi/2','pi','3\pi/2','pi'});
set(axes4,'YTick',[0 1 2 3 4],'YTickLabel',{'0','1','2','3','4'});
xlabel('$\theta$ (rad)','Interpreter','latex','FontSize',18);
ylabel('$\rho_{\rm NESS}(\theta)$(rad$^{-1})$', 'Interpreter','latex','FontSize',18);


set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',18,'Tickdir','in','TickLength',[0.02, 0.01],...
    'XTick',[0 pi  2*pi ],'XTickLabel',{'$0$','$\pi$','$2\pi$'})

%%BOX 3 
axes('Position',[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  
% 
% axes5 = axes('Parent',figure5);
% hold(axes5,'on');
ylim([-0.1 0.3]);
plot(tt(1:floor(L/kk)),C(1:floor(L/kk)),'LineWidth',3,'Color',col3)
hold on
plot([7.52,7.52],[-0.069,0.09],'--k', 'LineWidth',2)
plot([26.63,26.63],[-0.016,0.09],'--k', 'LineWidth',2)


yarr=0.27;
xrr = [0.2 0.366];    % adjust length and location of arrow 
yrr = [yarr yarr];      % adjust hieght and width of arrow
annotation('textarrow',xrr,yrr,'String',' ','FontSize',13,'Linewidth',2)

xrr = [0.2,0.163];    % adjust length and location of arrow 
yrr = [yarr yarr];      % adjust hieght and width of arrow
annotation('textarrow',xrr,yrr,'String',' ','FontSize',13,'Linewidth',2)

 text(15,0.12,'$\approx j^{-1}$','Interpreter','latex','FontSize',18)
% text(0.366,yarr,'A','FontSize',14)

% box(axes5,'on');
% set(axes5,'XTick',[0 10 20 30],'XTickLabel',{'0','10','20','30'});
% set(axes5,'YTick',[-0.1 0 0.1 0.2 0.3],'YTickLabel',{'-0.1','0','0.1','0.2','0.3'});
xlabel('$t$(s)','Interpreter','latex','FontSize',30);
ylabel('$C(t)$','Interpreter','latex','FontSize',30);

%  set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'XMinorTick','off','Tickdir','in','TickLength',[0.02, 0.01],...
%  'XTick',[0 10 20 30],'XTickLabel',{'$0$','$10$','$20$','$30$'})

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'XMinorTick','off','Tickdir','in','TickLength',[0.02, 0.01],...
  'XTick',[0 10 20 30],'XTickLabel',{'$0$','$10$','$20$','$30$'},'FontSize',25, 'xlim',[0 33.33])

%%BOX 4 
axes('Position',[bx1+xwi+bxm  0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  
% 


% axes6 = axes('Parent',figure6);
% hold(axes6,'on');
xlim([0 33.33]);
ylim([0 0.35]);
plot(tt(1:floor(L/kk)),C(1)-C(1:floor(L/kk)),'LineStyle','-.','LineWidth',3,'Color',[0 0 1])
hold on 

plot(tt(1:floor(L/kk)),B(1:floor(L/kk)),'Linestyle','--','LineWidth',3,'Color',[1 0 0])
plot(tt(1:floor(L/kk)),C(1)-C(1:floor(L/kk))-B(1:floor(L/kk)),'Linestyle','-','LineWidth',3,'Color',[0 0.498039215803146 0])
box(axes6,'on');
set(axes6,'XTick',[0 10 20 30],'XTickLabel',{'0','10','20','30'});
set(axes6,'YTick',[0 0.1 0.2 0.3 0.4],'YTickLabel',{'0','0.1','0.2','0.3','0.4'});
xlabel('$t$(s)','Interpreter','latex','FontSize',30);
ylabel('$\rho_{\rm NESS}(\theta)$ (rad$^{-1})$', 'Interpreter','latex','FontSize',30);
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'XMinorTick','off','Tickdir','in','TickLength',[0.02, 0.01],...
  'XTick',[0 10 20 30],'XTickLabel',{'$0$','$10$','$20$','$30$'},'FontSize',25, 'xlim',[0 33.33],'ylim',[0 0.35])


%%BOX 4 inset
axes('Position',[bx1+xwi+330 0 xwi/2 0]/Xpix + [0 by1+42 0 ywi/2.6]/Ypix);

axes7 = axes('Parent',figure7);
hold(axes7,'on');

plot(tt(1:floor(L/kk)),gamma*a*a*R1(1:floor(L/kk)),'LineStyle','-.','LineWidth',3,'Color',[0 0 1])
hold on
plot(tt(1:floor(L/kk)),gamma*a*a*R(1:floor(L/kk)),'Linestyle','-','LineWidth',3,'Color',[0 0.498039215803146 0])
box(axes7,'on');
set(axes7,'XTick',[0 10 20 30],'XTickLabel',{'0','10','20','30'});
set(axes7,'YTick',[0 2 4 6 8],'YTickLabel',{'0','2','4','6','8'});
xlabel('$t$(s)','Interpreter','latex','FontSize',30);
ylabel('$\gamma a^2 R(t)(\rm rad^{-2})$','Interpreter','latex','FontSize',30);

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'XMinorTick','off','Tickdir','in','TickLength',[0.02, 0.01],...
  'FontSize',18, 'xlim',[0 33.33],'ylim',[-1 8])


axes('Position',[(0) 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
hold on

xt = [120,270+xwi,120,270+xwi];
yt = [ Ypix-15,Ypix-15,by1+ywi+bym/5,by1+ywi+bym/5];

%  str = {'\bf a','\bf b','\bf c','\bf d'};
 str = {'\bf c','\bf d','\bf e','\bf f'};
text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])