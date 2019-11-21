%Compute storage G'(f) and loss G''(f) modulus of complex fluid 
%by Passive Microrheology from equilibrium Power Spectral Density (PSD)
%of embedded particle trapped by optical tweezers,
%using Kramers Kronig relation to find real part of response function
%from imaginary part, using Fluctuation-Dissipation theorem.
%A file containing the coordinates of the trapped article and the
%corresponding time is needed, and must be specified in the parameters 
%Filepath and Filename. This line could be replaced by load() is the data
%are provided as MATLAB files

% File CPyCl5mM.txt is needed to run PassiveMicrorheologyKramersKronig and get the plots of Figs. 2(c) and (d).
% This file contains the particle position (first column) and the corresponding time (second column).


close all
NFFT = 2^17; %number of points to compute PSD

%read files
Filepath = 'PassiveMicrorheologyData/';
Filename = 'CPyCl5mM';
Extension = '.txt';

T0=22; %bath temperature
kB=1.38e-23; %Boltzman constant
T0 =T0 + 273.16;
r = (2.0e-6)/2; %particle radius
filname = [Filepath Filename Extension];

%Read data
Data = dlmread(filname,'',1,0);


x = Data(:,1); %position in meters
t = Data(:,2); %time in seconds
fs = 1/(t(2)-t(1)); %sampling frequency


bx1 = 120;     % bordo a sinistra
xwi = 400;    % larghezza riquadro con funzione
bx2 = 30;     % bordino a destra

Xpix = 3*bx1+3*xwi+2*bx2;  % larghezza figura in pixel

by1 = 110;     % bordo in basso
ywi = 300;    % altezza riquadro con funzione
by2 = 50;     % bordo in alto

Ypix = by1+ywi+by2;  % larghezza figura in pixel



figure (1)
plot(t,x)
xlabel('t [seconds]')
ylabel('x [meters]')
title('Stochastic trajectory of trapped particle')

%Determine trapstifness by equipartition theorem
k=kB*T0/var(x) 

%Compute PSD
[f,spx] = spectave(detrend(x),fs,NFFT);

figure(2)
loglog(f,spx,'b')
hold on
xlabel('f [Hz]')
ylabel('PSD [m^2 / Hz]')
title('Equilibrium PSD of trapped particle')

%Smooth PSD using polynomial fit of degree D
D = 8;
logsx=log(spx(2:end));
logf = log(f(2:end));
coeffx=polyfit(logf,logsx,D);
spx1=exp(polyval(coeffx,(logf),D));
spx11=spx;
spx(2:end)=spx1;
spx12=spx;
figure(2)
loglog(f,spx12,'k')


%Compute imaginary (chi2) and real (chi1) parts of linear response function of x
%to perturbative force
chi2=pi*f.*spx/2/kB/T0; %Fluctuation-dissi�tion Theorem
chi11=KramersKronig(2*pi*f,chi2');
chi1=chi11';

%Compute storage (G1) and loss (G2) modulus from linear response function
G1 = chi1./(chi1.^2+chi2.^2)/6/pi/r - k/6/pi/r;
G2 = chi2./(chi1.^2+chi2.^2)/6/pi/r;


figure(3)
loglog(f,G1,'r o')
hold on
loglog(f,G2,'b s')
xlabel('f [Hz]')
ylabel('G�, G�� [Pa]')
title('Storage (circles) and loss (squares) modulus')




figure('Position',[10 20 Xpix Ypix]); % crea la figura
axes('Position',[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
loglog(f,spx11,'b')
hold on
xlabel('f [Hz]')
ylabel('PSD [m^2 / Hz]')
title('Equilibrium PSD of trapped particle')

%Smooth PSD using polynomial fit of degree D

loglog(f,spx12,'k')