clear all;
close all;


lwidth = 3;

filenm = ['data/1p7mbar-2p2mbar_10V_49kHz'];
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



figure('Position',[10 20 700 700]);

subplot(2,1,1);
%plot(time,real(exp(-1i*0.2664)*(m(:,2)+1i*m(:,3)))/calib,'o','MarkerSize',2);  %0.203
%hold all;
plot(time,imag(exp(-1i*0.2664)*(m(:,2)+1i*m(:,3)))/calib,'o','MarkerSize',2);
%xlim([940 1080]);
xlim([0 470]);
ylim([-9.5 5.5]*1E-8);
xlabel('$time (s)$', 'Interpreter','latex','FontSize',30);
ylabel('$amplitude (m_p)$','Interpreter','latex','FontSize',30);
title('switched off plasma after 8000 s', 'Interpreter','latex','FontSize',30);
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25)

subplot(2,1,2);
plot(time,imag(exp(-1i*0.2664)*(m(:,2)+1i*m(:,3)))/calib,'o','MarkerSize',2);
xlim([7950 9050]);
ylim([-5.5 5.5]*1E-8);
xlabel('$time (s)$','Interpreter','latex','FontSize',30);
ylabel('$amplitude (m_p)$', 'Interpreter','latex','FontSize',30);

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25)
%set(gcf, 'Position', [100, 100, 2049, 895]);
%set(gcf,'PaperPositionMode','auto');
%print(gcf,'-depsc',['CapSwitchOff_quadrature_meters','.eps']);


return


