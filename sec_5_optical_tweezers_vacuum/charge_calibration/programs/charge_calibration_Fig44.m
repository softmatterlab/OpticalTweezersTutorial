clear all;
close all;

col3=[242/255,177/255,223/255,0.3];
col4=[0.00,0.45,0.74,0.3];



filenm = ['../data/1p7mbar-2p2mbar_10V_49kHz'];

m = (csvread([filenm,'.csv']));    % read in data

time = m(:,1)-m(1,1);

% calculate amplitude of particle motion
% the DAC calibration gives 3270 bit/V
% including the amplifier in the filter/FB box we have 69000 bits/V
calib_DAC = 69000;
% the calibration of the particle motion at RT gives 1.4E8 bits/m
calib_RT = 1.4E8;

% calibration V/m is thus, minus sign comes from separate experiment,
% relevant to get sign of charge right
calib = -calib_RT/calib_DAC/sqrt(2)




xwi = 780;
xw2=300;
bx1 = 120;     % extra space at the left
bx2 = 100;     % extra space at the right
bxm = 160;     % extra space at the midle


Xpix=1400;   % total width
ywi = 200;    % length frame with function

by1 = 100;     % extra space below
by2 = 30;     % extra space up
bym = 130      %extra space at the midle
yw2=2*ywi+bym;
Ypix = by1+2*ywi+2*by2+bym;  % width in pixel


figure('Position',[100 200 Xpix Ypix]);

axes('Position',[bx1+xw2+bxm 0 xwi 0]/Xpix + [0 by1+ywi+bym 0 ywi]/Ypix);

%plot(time,real(exp(-1i*0.2664)*(m(:,2)+1i*m(:,3)))/calib,'o','MarkerSize',2);  %0.203
%hold all;
rectangle('Position',[-60,0,600,100],'FaceColor',col3)
hold on

rectangle('Position',[-60,-100,600,100],'FaceColor',col4)

scatter(time-58.63,imag(exp(-1i*0.2664)*(m(:,2)+1i*m(:,3)))/calib*1e9,5,'o','filled','MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth', 1);

plot([0 0],[-100 50],'--','LineWidth', 1.5, 'color','k')

annotation('textarrow',[0.52,0.49],[0.7,0.7],'String','HV turned on','Interpreter','latex','FontSize',22)
annotation('textarrow',[0.535,0.55],[0.85,0.81],'String','+1e','Interpreter','latex','FontSize',22)
annotation('textarrow',[0.58,0.595],[0.86,0.83],'String','+2e','Interpreter','latex','FontSize',22)
annotation('textarrow',[0.84,0.84],[0.70,0.73],'String','-1e','Interpreter','latex','FontSize',22)
annotation('textarrow',[0.89,0.89],[0.67,0.715],'String','-2e','Interpreter','latex','FontSize',22)


annotation('ellipse',[0.437 0.82 0.02 0.04 ],'LineWidth',1.5)
text(-39.5, 49,'+','Interpreter','latex','FontSize',25)
annotation('ellipse',[0.437 0.68 0.02 0.04 ],'LineWidth',1.5)
text(-37.1, -29,'_','Interpreter','latex','FontSize',25)


%xlim([940 1080]);
xlim([-58.63 400]);
ylim([-10 10]*1E1);
xlabel('$t$ (s)', 'Interpreter','latex','FontSize',30);
% ylabel('amplitude (nm)','Interpreter','latex','FontSize',30);
% title('switched off plasma after 8000 s', 'Interpreter','latex','FontSize',30);
box on
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25)

%%
axes('Position',[bx1+xw2+bxm 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);


rectangle('Position',[-50,0,1500,50],'FaceColor',col3)
hold on

rectangle('Position',[-50,-50,1500,50],'FaceColor',col4)


scatter(time-7974,imag(exp(-1i*0.2664)*(m(:,2)+1i*m(:,3)))/calib*1e9,5,'o','filled','MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth', 1);
plot([0 0],[-50 50],'--','LineWidth', 1.5, 'color','k')

annotation('textarrow',[0.46,0.43],[0.2,0.2],'String','HV turned off','Interpreter','latex','FontSize',22)



xlim([7950-7974 1000]);
ylim([-5 5]*1E1);
xlabel('$t$ (s)','Interpreter','latex','FontSize',30);
% ylabel('amplitude (nm)', 'Interpreter','latex','FontSize',30);
box on

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25)
%set(gcf, 'Position', [100, 100, 2049, 895]);
%set(gcf,'PaperPositionMode','auto');
%print(gcf,'-depsc',['CapSwitchOff_quadrature_meters','.eps']);

fsample = 625000;
% fsize = 15;

PSDall = zeros(32769,3);

fprintf('#################\n###########\n remember to enter correct pressure\n ##########\n###############\n');
% gas pressure during measurement in Pascal, which is mbar x 100
Pgas = 1.9E2;

% get calibration to convert voltage to particle position

% first, get all data from calibration measurements
for j=0:19
    if j<10
        jstring = ['0',num2str(j)];
    else
        jstring = num2str(j);
    end;
    %cd 'TimeTraces'
    filenm = ['../data/17233_p=1p9mbar_V=10V+0.1V_f=45p3kHz_0',jstring];
    %M = csvread([filenm,num2str(j),'.tt']);
    m = (csvread([filenm,'.tt']))';    % read in data
    %cd ..
    [PSDm, absFFTout, FFTm, fm]=FFTnormalized(m, fsample);
    PSDall = PSDall+PSDm;
    var_m(j+1,:) = var(m);
end;

varMean = mean(var_m)   % this variance can be used to gauge amplitudes of signal^2 over time
PSDmean = PSDall/(j+1);



kb = 1.38E-23;  % Boltzmann constant
T = 293;    % room temperature in K
rhoGlass = 2.2E3; % glass density in kg/m^3
eta = 18.27E-6;
dm = 0.372E-9;


% rpart = 132E-9/2;   % particle radius
% Vpart = 4*pi/3*rpart^3;
% mpart = rhoGlass*Vpart;


% fit Lorentzian to data
% select range of data to fit to
fitRangeHz2 = [40,50]*1E3;   % fit range in Hertz
fitRange = round(fitRangeHz2/(fm(2)-fm(1)));

fitRange_excludeHz = [45200,45400];
fitRange_ex = round(fitRange_excludeHz/(fm(2)-fm(1)));

dataToFit = PSDmean(:,1);
dataToFit(fitRange_ex(1)+1:fitRange_ex(2)) = dataToFit(fitRange_ex(1));
dataToFit = dataToFit(fitRange(1)+1:fitRange(2),1);

startingGuess = [1E14 , 1000, 45000, 0];
%[estimates, sse, model] = fitLorentz(fm', PSDmean(:,2), startingGuess);
[estimates, sse, model] = fitLorentz((fitRange(1)+1:fitRange(2))'*(fm(2)-fm(1)), dataToFit, startingGuess);
estimates1 = estimates;
[sse, FittedCurve] = model(estimates);

%calib1 = sqrt(kb*T/(pi*estimates(1)*mpart))
freq1 = estimates(3)

% derive the mass from the particle radius using Erik's instructions
% particle radius
r1 = 0.619*9*pi/sqrt(2)*eta*dm^2*Pgas/(rhoGlass*kb*T*2*pi*abs(estimates(2)))
% get derived particle mass
m1 = 4*pi/3*r1^3*rhoGlass;
calibErik1 = 2*pi*sqrt(abs(estimates(1))*pi*m1/(2*abs(estimates(2))*kb*T))
% get calibration using equipartition
calibEquip1 = sqrt(m1*(2*pi*freq1)^2*varMean(1)/kb/T)
%%

% axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
axes('Position',[bx1+xw2+110 0 xwi 0]/Xpix + [0 by1+bym+30 0 ywi]/Ypix,'visible','off');
ylabel('quadrature amplitude (nm)', 'Interpreter','latex','FontSize',30,'visible','on');

% plot for paper
%%
axes('Position',[bx1 0 xw2 0]/Xpix + [0 by1 0 yw2]/Ypix);
fm_short = reshape(fm(1:end-1),4,[]);
fm_short = fm_short(1,:);
PSDmean_short = reshape(PSDmean(1:end-1,1),4,[]);
PSDmean_short = PSDmean_short(1,:);
%p=semilogy(fm,PSDmean(:,1)/calibErik1^2,(fitRange(1)+1:fitRange(2))'*(fm(2)-fm(1)),FittedCurve/calibErik1^2);
p=semilogy(fm*1e-3,PSDmean(:,1)/calibErik1^2*1e18,'o','color','k','LineWidth',1)
hold on
p1=semilogy((fitRange(1)+1:fitRange(2))'*(fm(2)-fm(1))*1e-3,FittedCurve/calibErik1^2*1e18,'LineWidth',2,'color','r');
% p=semilogy(fm,PSDmean(:,1)/calibErik1^2,'o',(fitRange(1)+1:fitRange(2))'*(fm(2)-fm(1)),FittedCurve/calibErik1^2,'LineWidth',2,'color','k');
%p=semilogy(fm_short,PSDmean_short/calibErik1^2,'o');

plot([4.53e1 4.53e1],[1e-2 1e+8],'--','LineWidth', 1.5, 'color','k')

%  text(4.3e4,2e-16,'f$_d$','Interpreter','latex','FontSize',22,'BackgroundColor', 'w')
annotation('textarrow',[0.19,0.22],[0.7,0.68],'String','$f_{\rm d}$ ','Interpreter','latex','FontSize',22)
xlim([40 48])
ylim([2E-1 3E+3]);
% set(p,'LineWidth',3);
% set(gca,'FontSize',fsize);
xlabel('$f$ (kHz)','Interpreter','latex','FontSize',30);
ylabel('PSD (nm$^2$  Hz$^{-1})$','Interpreter','latex','FontSize',30);
% title(['axis 1, Lorentzian fit, d1=',num2str(2*round(r1*1E9)),'nm']);
% set(gcf,'PaperPositionMode','auto');
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'Xtick',[40 42 44 46 48 ])


%%
axes('Position',[0 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix);


xt = [bx1,bx1+xw2+bxm,bx1+xw2+bxm];
yt = [by1+yw2+by2,by1+yw2+by2,by1+ywi+by2];

str = {'\bf a','\bf b','\bf c'};


text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])


%saveas(gcf,'Fig43.eps','epsc')
%saveas(gcf,'Fig43.fig')




