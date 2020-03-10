%This program computes the work (Work), the heat (Heat), the
%variation of the potential energy (Energy) and the total entropy production (totEntropy), 
%measured over a time interval of duration tau, for a spherical particle (radius r) driven 
%in water (viscosity gamma, temperature T0) by an optical trap (stiffness k) moving in  
%a straightline at constant velocity v. The stochastic thermodynamics quantities are determined  
%from the time (t) evolution of the particle position (x) and
%the time evolution of the position of the optical trap (xtrap). 
%To this end, first the numerical values of k and v are determined from x and xtrap,
%respectively, and then they are used as input parameters within the formalism of Stochastic
%Thermodynamics. The probability density functions of the stochastic
%thermodynamic quantities are also computed, as well as their mean values
%and standars deviations. Finally, the integral fluctuation theorem is
%checked for the total entropy production at different measurement times
%tau.

clear all
close all

%read data
Filepath = 'Fig29/';
Filename = 'StochThermoTrajectory_data';
Extension = '.txt';

T0=22; %bath temperature
kB=1.38e-23; %Boltzman constant
T0 =T0 + 273.16;
r = (2.73e-6)/2; %particle radius (m)
eta = 0.95e-3; %water viscosity at T0
gamma = 6*pi*eta*r; %friction coefficient


filname = [Filepath Filename Extension]; 

%Read data
Data = dlmread(filname,'',1,0);

xtrap = Data(:,1); %trap position (microns)
x = Data(:,2); %particle position (microns)
t = Data(:,3); %time (seconds)

fs = 1/(t(2)-t(1)); %sampling frequency

figure(1)
figure1 = figure(1);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
xlim(axes1,[25 30]);
ylim(axes1,[-1.5 1.5]);
plot(t,x)
hold on
plot(t,xtrap,'Linestyle','--')
box(axes1,'on');
set(axes1,'XTick',[]);
set(axes1,'YTick',[-1 0 1],'YTickLabel',{'-1','0','1'});
ylabel({'{\itx} [\mum]'});
xp=x;
xtrapp=xtrap;

dx = xtrap - x; %-(particle position relative to trap position) 
dxg=dx;
figure(2)
figure2 = figure(2);
axes2 = axes('Parent',figure2);
hold(axes2,'on');
xlim(axes2,[25 30]);
ylim(axes2,[-0.25 0.25]);
plot(t,dx)
box(axes2,'on');
set(axes2,'XTick',[25 26 27 28 29 30],'YTickLabel',{'25','26','27','28','29','30'});
set(axes2,'YTick',[-0.2 0 0.2],'YTickLabel',{'-0.2','0','0.2'});
xlabel({'{\itt}[s]'});
ylabel({'\lambda({\itt}) - {\itx} [\mum]'});

%Determine trap stiffness k and driving velocity v, first convert length-scale units into
%meter
x = 1e-6*x;
xtrap = 1e-6*xtrap;
dx = 1e-6*dx;

k = kB*T0/var(dx) %trap stiffness (N/m)
coeff = polyfit(t,xtrap,1);
v = coeff(1) %driving velocity imposed by the moving trap(m/s)

%Compute fluctuations (Work) of the work (w) done by the trap, the fluctuations (Heat) of the heat (q) 
%dissipated into the fluid, the fluctuations (Energy) of the variations of the potentiañ energy (du),
%and the fluctuations (totEntropy) of total entropy for different measurement times tau. 

%-Determine the probability densities pdfw, pdfq and pdfdu of w, q and du,
 %respectively, as well as the mean values wm, qm, dum and the standard
%deviations ws, qs, dus.

%-Compute the average IFT = <exp(-totEntropy/kB)> involved in the Integral
 %Fluctuation Theorem for different observation times tau

t0 = 0; %initial point
t1 = length(dx); %lasti point
Ntau = 6; %number of different values of tau to compute pdfs 
for tt = 1:1:Ntau
    ntau = 5*tt; %five times the temporal resolution
    tau(tt) = ntau/fs; %(observation/measurement) time interval
    for ii=(t0:(t1 - ntau - 2))
        Work(ii-t0+1) = k*v*(sum(dx(ii+1:ii+ntau+1))-0.5*(dx(ii+1)+dx(ii+ntau+1)))/fs;   
        Energy(ii-t0+1) = 0.5*k*(dx(ii+ntau+1)^2-dx(ii+1)^2);
        int = -(dx(ii+1:ii+1+ntau)+dx(ii+2:ii+2+ntau)).*(x(ii+2:ii+2+ntau)-x(ii+1:ii+1+ntau));
        DissHeat = -0.5*k*(sum(int)-0.5*(int(1)+int(end))); 
        Heat(ii-t0+1) = DissHeat;
        totEntropy(ii-t0+1) = DissHeat/T0 + 0.5*k/T0*((-dx(ii+ntau+1)+gamma*v/k)^2-(-dx(ii+1)+gamma*v/k)^2);
        clear int DissHeat
    end
    [PDFW,W]=histnorm(Work, 100, 0);
    [PDFdU,dU]=histnorm(Energy, 100, 0);
    [PDFQ,Q]=histnorm(Heat, 100, 0);
    wm(tt)=mean(Work);
    ws(tt)=std(Work);
    dum(tt)=mean(Energy);
    dus(tt)=std(Energy);
    qm(tt)=mean(Heat);
    qs(tt)=std(Heat);
    w(:,tt)=W;
    pdfw(:,tt)=PDFW;
    du(:,tt)=dU;
    pdfdu(:,tt)=PDFdU;
    q(:,tt)=Q;
    pdfq(:,tt)=PDFQ;
    %Verify Integral Fluctuation Theorem for total entropy production
    IFT(tt) = mean(exp(-totEntropy/kB)); 
    clear W PDFW dU PDFdU Q PDFQ Work Energy Heat totEntropy
end

%Add zeros to the mean and standard deviations for the trivial case tau = 0:
tau = [0 tau];
wm = [0 wm];
ws = [0 ws];
qm = [0 qm];
qs = [0 qs];
dum = [0 dum];
dus = [0 dus];
IFT = [exp(-0) IFT];

figure(3)
figure3 = figure(3);
axes3 = axes('Parent',figure3);
hold(axes3,'on');
xlim(axes3,[0 0.21]);
ylim(axes3,[-0.001 0.4]);
plot(tau,wm/kB/T0,'LineStyle','none','MarkerFaceColor',[0 0 1],'Marker','o','Color',[0 0 1],'MarkerSize',10)
hold on
plot(tau,qm/kB/T0,'LineStyle','none','MarkerFaceColor',[1 0 0],'MarkerSize',8,'Marker','square',...
    'Color',[1 0 0])
plot(tau,dum/kB/T0,'LineStyle','none','MarkerFaceColor',[0 0.498039215803146 0],'Marker','diamond',...
    'Color',[0 0.498039215803146 0],'MarkerSize',10)
plot(tau,gamma*v^2*tau/kB/T0,'LineWidth',1,'Color',[0 0 0],'LineStyle','-')
box(axes3,'on');
set(axes3,'XTick',[0 0.05 0.1 0.15 0.2],'YTickLabel',{'0','0.05','0.1','0.15','0.2'});
set(axes3,'YTick',[0 0.1 0.2 0.3 0.4],'YTickLabel',{'0','0.1','0.2','0.3','0.4'});
xlabel({'\tau [s]'});
ylabel({'mean [{\itk_BT}]'});

ttt=0:0.01:0.21;
figure(4)
figure4 = figure(4);
axes4 = axes('Parent',figure4);
hold(axes4,'on');
xlim(axes4,[0 0.21]);
ylim(axes4,[0 1.5]);
plot(tau,ws/kB/T0,'LineStyle','none','MarkerFaceColor',[0 0 1],'Marker','o','Color',[0 0 1],'MarkerSize',10)
hold on
plot(tau,qs/kB/T0,'LineStyle','none','MarkerFaceColor',[1 0 0],'MarkerSize',8,'Marker','square',...
    'Color',[1 0 0])
plot(tau,dus/kB/T0,'LineStyle','none','MarkerFaceColor',[0 0.498039215803146 0],'Marker','diamond',...
    'Color',[0 0.498039215803146 0],'MarkerSize',10)
plot(ttt,sqrt(2*gamma*v^2/kB/T0)*sqrt(ttt + gamma/k*exp(-k*ttt/gamma)-gamma/k),'LineWidth',1,'Color',[0 0 0],'LineStyle','-')
box(axes4,'on');
set(axes4,'XTick',[0 0.05 0.1 0.15 0.2],'YTickLabel',{'0','0.05','0.1','0.15','0.2'});
set(axes4,'YTick',[0 0.5 1 1.5],'YTickLabel',{'0','0.5','1','1.5'});
xlabel({'\tau [s]'});
ylabel({'standard deviation [{\itk_BT}]'});

%Determine Gaussian fits of work probability density function
ww = (-3:0.01:3)*kB*T0;
for tt = 1:Ntau
   pdfGauss(:,tt) = exp(-0.5*(ww - wm(tt)).^2/ws(tt)^2)/sqrt(2*pi)/ws(tt);
end

figure(5)
figure5 = figure(5);
axes5 = axes('Parent',figure5);
hold(axes5,'on');
xlim(axes5,[-2.5 3]);
ylim(axes5,[0 2]);
plot(w/kB/T0,pdfw*kB*T0,'Marker','.','Linestyle','none','MarkerSize',9)
hold on
plot(ww/kB/T0,pdfGauss*kB*T0,'k')
box(axes5,'on');
set(axes5,'XTick',[-2 -1 0 1 2 3],'YTickLabel',{'-2','-1','0','1','2','3'});
set(axes5,'YTick',[0 0.5 1 1.5 2],'YTickLabel',{'0','0.5','1','1.5','2'});
xlabel({'{\itW}_{\tau} [{\itk_BT}]'});
ylabel({'{\itP}({\itW}_{\tau}) [({\itk_BT})^{-1}]'});


figure(6)
figure6 = figure(6);
axes6 = axes('Parent',figure6);
hold(axes6,'on');
xlim(axes6,[-10 15]);
ylim(axes6,[0 1]);
plot(q/kB/T0,pdfq*kB*T0,'LineWidth',1)
box(axes6,'on');
set(axes6,'XTick',[-10 -5 0 5 10 15],'YTickLabel',{'-10','-5','0','5','10','15'});
set(axes6,'YTick',[0 0.2 0.4 0.6 0.8 1],'YTickLabel',{'0','0.2','0.4','0.6','0.8','1'});
xlabel({'{\itQ}_{\tau} [{\itk_BT}]'});
ylabel({'{\itP}({\itQ}_{\tau}) [({\itk_BT})^{-1}]'});

figure(7)
figure7 = figure(7);
axes7 = axes('Parent',figure7);
hold(axes7,'on');
xlim(axes7,[0 0.21]);
ylim(axes7,[0 2]);
plot(tau,IFT,'MarkerFaceColor',[0 0 0],'MarkerSize',10,'Marker','pentagram',...
    'LineStyle','none',...
    'Color',[0 0 0])
box(axes7,'on');
set(axes7,'XTick',[0 0.1 0.2],'YTickLabel',{'0','0.1','0.2'});
set(axes7,'YTick',[0 1 2],'YTickLabel',{'0','1','2'});
xlabel({'\tau [s]'});
ylabel({'\langle exp(-\Delta{\itS}_{\tau}/{\itk_B}) \rangle'});


figure(8)
figure8 = figure(8);
axes8 = axes('Parent',figure8);
hold(axes8,'on');
xlim(axes8,[-10 10]);
ylim(axes8,[0 1.2]);
plot(du/kB/T0,pdfdu*kB*T0,'LineWidth',1)
box(axes8,'on');
set(axes8,'XTick',[-10 -5 0 5 10],'YTickLabel',{'-10','-5','0','5','10'});
set(axes8,'YTick',[0 0.2 0.4 0.6 0.8 1 1.2],'YTickLabel',{'0','0.2','0.4','0.6','0.8','1','1.2'});
xlabel({'\Delta{\itU}_{\tau} [{\itk_BT}]'});
ylabel({'{\itP}(\Delta{\itU}_{\tau}) [({\itk_BT})^{-1}]'});





%% 
%% 

xwi = 315;    % width of the plot square
bx1 = 150;     % extra space at the left
bx2 = 20;     % extra space at the right
bxm = 140;     % extra space at the midle


Xpix=1400;   % total width 
ywi = 170;    % length frame with function

by1 = 50;     % extra space below
by2 = 30;     % extra space up
bym = 80      %extra space at the midle
Ypix = by1+2*ywi+by2+bym;  % width in pixel


col3=[0.00,0.45,0.74];

figure('Position',[10 20 Xpix 2*Ypix]); % generate te figure



%%BOX 1
axes('Position',[bx1 0 xwi 0]/Xpix + [0 by1+ywi+(1-0.65)*ywi+bym 0 0.65*ywi]/Ypix);  

hold on;
xlim([25 30]);
ylim([-1.5*1e3 1.5*1e3]);
plot(t,xp*1e3)
hold on
plot(t,xtrapp*1e3,'Linestyle','--')
box on;
set(axes1,'XTick',[]);
set(axes1,'YTick',[-1 0 1],'YTickLabel',{'-1','0','1'});
ylabel('$x$ (nm)','Interpreter','latex','FontSize',30);


set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'YMinorTick','off','Tickdir','in','XTickLabel',{[]})

%-(particle position relative to trap position) 
% 'TickLength',[0.02, 0.01],...
%         'YTick',[0 0.25 0.50 0.75 1 ]
    


axes('Position',[bx1 0 xwi 0]/Xpix + [0 by1+ywi+bym 0 0.3*ywi]/Ypix);  
xlim([25 30]);
ylim([-0.25*1e3 0.25*1e3]);
plot(t,dxg*1e3)
box on;
%  set('XTick',[25 26 27 28 29 30],'YTickLabel',{'25','26','27','28','29','30'});
% set('YTick',[-0.2 0 0.2],'YTickLabel',{'-0.2','0','0.2'});
xlabel('$t$ (s)','Interpreter','latex','FontSize',30);
ylabel('$\lambda(t) - x$ (nm)','Interpreter','latex','FontSize',30);
xlim([25 30]);
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'YMinorTick','off','Tickdir','in',...
    'XTick',[25 26 27 28 29 30])


%%BOX 2
axes('Position',[bx1+xwi+bxm 0 xwi 0]/Xpix + [0 by1+ywi+bym 0 ywi]/Ypix);  

xlim([0 0.21]);
ylim([-0.001 0.4]);
plot(tau,wm/kB/T0,'LineStyle','none','MarkerFaceColor',col3,'Marker','o','Color',col3,'MarkerSize',10)
hold on
plot(tau,qm/kB/T0,'LineStyle','none','MarkerFaceColor',[1 0 0],'MarkerSize',8,'Marker','square',...
    'Color',[1 0 0])
plot(tau,dum/kB/T0,'LineStyle','none','MarkerFaceColor',[0 0.498039215803146 0],'Marker','diamond',...
    'Color',[0 0.498039215803146 0],'MarkerSize',10)
plot(tau,gamma*v^2*tau/kB/T0,'LineWidth',1,'Color',[0 0 0],'LineStyle','-')
box(axes3,'on');
set(axes3,'XTick',[0 0.05 0.1 0.15 0.2],'YTickLabel',{'0','0.05','0.1','0.15','0.2'});
set(axes3,'YTick',[0 0.1 0.2 0.3 0.4],'YTickLabel',{'0','0.1','0.2','0.3','0.4'});
xlabel('$\tau$ (s)','Interpreter','latex','FontSize',30);
ylabel('mean $(k_{\rm B}T)$','Interpreter','latex','FontSize',30);

ylabel('$\langle W_{\tau} \rangle\ (k_{\rm B}T) $','Interpreter','latex','FontSize',30);

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'YMinorTick','off','Tickdir','in')



%%BOX 3
axes('Position',[bx1+2*(xwi+bxm) 0 xwi 0]/Xpix + [0 by1+ywi+bym 0 ywi]/Ypix);  

xlim([0 0.21]);
ylim([0 1.5]);
plot(tau,ws/kB/T0,'LineStyle','none','MarkerFaceColor',col3,'Marker','o','Color',col3,'MarkerSize',10)
hold on
plot(tau,qs/kB/T0,'LineStyle','none','MarkerFaceColor',[1 0 0],'MarkerSize',8,'Marker','square',...
    'Color',[1 0 0])
plot(tau,dus/kB/T0,'LineStyle','none','MarkerFaceColor',[0 0.498039215803146 0],'Marker','diamond',...
    'Color',[0 0.498039215803146 0],'MarkerSize',10)
plot(ttt,sqrt(2*gamma*v^2/kB/T0)*sqrt(ttt + gamma/k*exp(-k*ttt/gamma)-gamma/k),'LineWidth',1,'Color',[0 0 0],'LineStyle','-')
box(axes4,'on');
set(axes4,'XTick',[0 0.05 0.1 0.15 0.2],'YTickLabel',{'0','0.05','0.1','0.15','0.2'});
set(axes4,'YTick',[0 0.5 1 1.5],'YTickLabel',{'0','0.5','1','1.5'});
xlabel('$\tau$ (s)','Interpreter','latex','FontSize',30);
ylabel('$\sigma(\langle W_{\tau} \rangle)\ (k_{\rm B}T)$','Interpreter','latex','FontSize',30);

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'YMinorTick','off','Tickdir','in')


%%BOX 4

axes('Position',[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  

axes5 = axes('Parent',figure5);
hold(axes5,'on');
xlim(axes5,[-2.5 3]);
ylim(axes5,[0 2]);
plot(w/kB/T0,pdfw*kB*T0,'Marker','.','Linestyle','none','MarkerSize',9)
hold on
plot(ww/kB/T0,pdfGauss*kB*T0,'k')
box(axes5,'on');
set(axes5,'XTick',[-2 -1 0 1 2 3],'YTickLabel',{'-2','-1','0','1','2','3'});
set(axes5,'YTick',[0 0.5 1 1.5 2],'YTickLabel',{'0','0.5','1','1.5','2'});
xlabel('$W_{\tau}\ (k_{\rm B}T)$','Interpreter','latex','FontSize',30) ;
ylabel('$P(W_{\tau}) (k_{\rm B}T)^{-1}$','Interpreter','latex','FontSize',30);

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'YMinorTick','off','Tickdir','in')


%%BOX 5
axes('Position',[bx1+xwi+bxm 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  

plot(q/kB/T0,pdfq*kB*T0,'LineWidth',1)
box on;
xlabel('$Q_{\tau}\ (k_{\rm B}T)$','Interpreter','latex','FontSize',30) ;
ylabel('$P(Q_{\tau})\  (k_{\rm B}T)^{-1}$','Interpreter','latex','FontSize',30);
xlim([-10 15]);
ylim([0 1]);
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'YMinorTick','off','Tickdir','in',...
    'XTick',[-10 -5 0 5 10 15],'YTickLabel',{'-10','-5','0','5','10','15'},'YTick',[0 0.2 0.4 0.6 0.8 1],'YTickLabel',...
    {'0','0.2','0.4','0.6','0.8','1'})


%%BOX 5 inset

axes('Position',[bx1+xwi+bxm+190 0 xwi/3 0]/Xpix + [0 by1+ywi/2 0 ywi/2.5]/Ypix);  


plot(tau,IFT,'MarkerFaceColor',[0 0 0],'MarkerSize',10,'Marker','pentagram',...
    'LineStyle','none',...
    'Color',[0 0 0])
box on;

ylabel('$\langle exp(-\Delta S_{\tau}/k_{\rm B}) \rangle $','Interpreter','latex','FontSize',30);
xlabel('$\tau$ (s)','Interpreter','latex','FontSize',30);
xlim([0 0.21]);
ylim([0 2]);
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',18,'YMinorTick','off','Tickdir','in',...
    'XTick',[0 0.1 0.2],'YTickLabel',{'0','0.1','0.2'},'YTick',[0 1 2],'YTickLabel',{'0','1','2'})

%%BOX 6
axes('Position',[bx1+2*(xwi+bxm) 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  

xlim([-10 10]);
ylim([0 1.2]);
plot(du/kB/T0,pdfdu*kB*T0,'LineWidth',1.5)
box on;
% set('XTick',[-10 -5 0 5 10],'YTickLabel',{'-10','-5','0','5','10'});
% set('YTick',[0 0.2 0.4 0.6 0.8 1 1.2],'YTickLabel',{'0','0.2','0.4','0.6','0.8','1','1.2'});
% xlabel({'\Delta{\itU}_{\tau} [{\itk_BT}]'});
% ylabel({'{\itP}(\Delta{\itU}_{\tau}) [({\itk_BT})^{-1}]'});

xlabel('$\Delta U_{\tau}\ (k_{\rm B}T)$','Interpreter','latex','FontSize',30) ;
ylabel('$P(\Delta U_{\tau}) (k_{\rm B}T)^{-1}$','Interpreter','latex','FontSize',30);

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'YMinorTick','off','Tickdir','in',...
    'XTick',[-10 -5 0 5 10],'YTickLabel',{'-10','-5','0','5','10'},...
    'YTick',[0 0.2 0.4 0.6 0.8 1 1.2],'YTickLabel',{'0','0.2','0.4','0.6','0.8','1','1.2'})



%% Indexes 
axes('Position',[(0) 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
hold on

xt = [bx1,bx1+xwi+bxm,bx1+2*(xwi+bxm),bx1,bx1+xwi+bxm,bx1+2*(xwi+bxm)];
yt = [by1+2*ywi+bym+15,by1+2*ywi+bym+15,by1+2*ywi+bym+15,by1+ywi+15,by1+ywi+15,by1+ywi+15];

 str = {'\bf a','\bf b','\bf c','\bf d','\bf e','\bf f'};


text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])