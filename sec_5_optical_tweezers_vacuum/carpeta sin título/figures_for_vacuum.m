%%FIGURES FOR RUBEN
clear all
close all

%FIG1
clear all
figure(1)
load('data1_vacuum.mat')

subplot(1,4,1);
loglog(pressure, gamma_tot/(2*pi), 'LineWidth', 2);
hold on
loglog(pressure, gamma_tot_2_bath/(2*pi),'--', 'LineWidth', 2)
xlabel('pressure (mBar)')
ylabel('damping rate $\gamma/(2\pi)$ (Hz)')
hold off


subplot(1,4,2);
loglog(pressure, Gamma/kB, 'LineWidth', 2);
hold on
loglog(pressure, Gamma_2_bath/kB,'--', 'LineWidth', 2)
xlabel('pressure (mBar)')
ylabel('heating rate $\Gamma/k_B$ (K/s)')
hold off


subplot(1,4,3);
loglog(pressure, T_FB, 'LineWidth', 2);
hold on
loglog(pressure,  T_FB_2_bath,'--', 'LineWidth', 2)
xlabel('pressure (mBar)')
ylabel('COM temperature (K)')
hold off

subplot(1,4,4);
loglog(pressure, T_int);
hold on
xlabel('pressure (mBar)')
ylabel('internal temperature (K)')
hold off

clear all
figure(2)

load('data2_vacuum_0.mat')
loglog(frequency/1000, PSD/g0, 'LineWidth', 2)
hold on 
loglog(frequency/1000, Lorentz,'--', 'LineWidth', 2, 'Color', 'black')
plot([x0 x0],[ymin ymax])

load('data2_vacuum_1.mat')
loglog(frequency/1000, PSD/g0, 'LineWidth', 2)
hold on 
loglog(frequency/1000, Lorentz,'--', 'LineWidth', 2, 'Color', 'black')
plot([x0 x0],[ymin ymax])

load('data2_vacuum_2.mat')
loglog(frequency/1000, PSD/g0, 'LineWidth', 2)
hold on 
loglog(frequency/1000, Lorentz,'--', 'LineWidth', 2, 'Color', 'black')
plot([x0 x0],[ymin ymax])
xlabel('$Frequency f (kHz)$', 'Interpreter', 'Latex')
ylabel('$Normalized PSD\n \hat{S}_{vv}(f)/g (bit^2/Hz^2)$', 'Interpreter', 'Latex')
hold off

clear all
figure(3)
load('data3_vacuum.mat')
subplot(2,1,1)
semilogy(f1, PSD1)
hold on
plot([xx0 xx0],[1e0 1e3])
plot([xx1 xx1],[1e0 1e3])
plot([xx2 xx2],[1e0 1e3])
ylabel('PSD')


subplot(2,1,2)
semilogy(f2, PSD2)
hold on
plot([XV25000 XV25000],[1e-2 1e3])
plot([XV45000 XV45000],[1e-2 1e3])
plot([XV50000 XV50000],[1e-2 1e3])
ylabel('PSD')

xlabel('frequency[Hz]')