%This program computes the real and the imaginary part of the complex
%modulus of a viscoelastic material, measured by Active Microrheology.
%The particle radius and trap stiffness must be specified as an input
%A file containing the coordinates of the trapped article, the sinusoidal 
%perturbative force at a single frequency and the corresponding time, 
%is needed, and must be specified in the parameters 
%filename. 



clear all
close all

xwi = 550;    % width of the plot square
bx1 = 100;     % extra space at the left
bx2 = 50;     % extra space at the right

Xpix = 2*xwi+2*bx1+2*bx2;  % total
Xpix=1400;
ywi = 450;    % length riquadro con funzione

by1 = 75;     % extra space below
by2 = 30;     % extra space up

Ypix = 1*by1+1*ywi+1*by2;  % larghezza figura in pixel
%number of bins of the histogram, if not set default is 50
figure('Position',[10 20 Xpix Ypix]);
R = 1.0e-6; %particle radius
k = 1e-6; %trap stiffness (previously measure ny equipartition in thermal equilibrium)

%Write driving frequencies at which sinusoidal motion of the optial trap
%was applied, they must be specified in the file name
FF = [0.05 0.1 0.2 0.4 0.8 1.0 1.6 3.2 6.4 12.8 25.6];

Filnam ='ActiveMicrorheologyData/5mMdrive_k1e-06_freq_'
positionintefig1=[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix;
axes( 'Position',positionintefig1);
for i = 1:length(FF)

ff = FF(i);
    
filename = [Filnam  num2str(ff) '.mat']; 

load(filename)

%read data
fs=data.fs; %sampling frequency
x=data.x; %particle position
fd=data.fd; %perturbative sinusoidal force
t=data.t; %time
if i==3
 xsel=x;
 tsel=t;
end
nfft = 2^17; %number of point to compute cross power spectral density

[pxx,f] = cpsd(x,x,[],[],nfft,fs);
[pxfd,f] = cpsd(fd,x,[],[],nfft,fs);

invchi = pxfd./pxx; %1/Fourier transform of response function of x to fd

reinvchi = real(invchi); %real part of 1/Fourier transform of response function of x to fd
iminvchi = imag(invchi); %imaginary part of 1/Fourier transform of response function of x to fd


loglog(f,pxx)
hold on

[mm,I]=(min(abs(ff - f))); 
freqd = f(I);

frqd(i) = freqd;
G1(i) = (reinvchi(I) - k)/6/pi/R %storage modulus
G2(i) = iminvchi(I)/6/pi/R %loss modulus

end

xlabel('$f [Hz]$', 'Interpreter','Latex', 'FontSize',30)
ylabel('$PSD [m^2/Hz]$', 'Interpreter','Latex', 'FontSize',30)
set(gca,'TickLabelInterpreter','Latex', 'linewidth',1.5,'FontSize',25);
hold off

positionintefig2=[bx1+50 0 xwi/2 0]/Xpix + [0 3*by1/2 0 ywi/3]/Ypix;
axes( 'Position',positionintefig2);
plot(tsel, xsel)
xlim([0 40])

%%
positionintefig3=[2*bx1+xwi+bx2 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix;
axes( 'Position',positionintefig3);
loglog(frqd,G1,'r o')
hold on
loglog(frqd,G2,'s b')
xlabel('f [Hz]', 'Interpreter','Latex', 'FontSize',30)
ylabel('G(f),G''(f) [Pa]', 'Interpreter','Latex', 'FontSize',30)
%title('Storage (circles) and loss (squares) modulus of material obtained bt Active Microrheology')
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01]);
%%


%we add the data of the passive case
NFFT = 2^17;
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


% figure(1)
% plot(t,x)
% xlabel('t [seconds]')
% ylabel('x [meters]')
% title('Stochastic trajectory of trapped particle')

%Determine trapstifness by equipartition theorem
k=kB*T0/var(x) 

%Compute PSD
[f,spx] = spectave(detrend(x),fs,NFFT);


%loglog(f,spx,'b')
%hold on

%title('Equilibrium PSD of trapped particle')

%Smooth PSD using polynomial fit of degree D
D = 8;
logsx=log(spx(2:end));
logf = log(f(2:end));
coeffx=polyfit(logf,logsx,D);
spx1=exp(polyval(coeffx,(logf),D));
spx(2:end)=spx1;


%loglog(f,spx,'k')


%Compute imaginary (chi2) and real (chi1) parts of linear response function of x
%to perturbative force
chi2=pi*f.*spx/2/kB/T0; %Fluctuation-dissi�tion Theorem
chi11=KramersKronig(2*pi*f,chi2');
chi1=chi11';

%Compute storage (G1) and loss (G2) modulus from linear response function
G1 = chi1./(chi1.^2+chi2.^2)/6/pi/r - k/6/pi/r;
G2 = chi2./(chi1.^2+chi2.^2)/6/pi/r;



loglog(f,G1,'r --')

loglog(f,G2,'b --')
%xlabel('f [Hz]')
%ylabel('G�, G�� [Pa]')
%title('Storage (circles) and loss (squares) modulus')

xlim([1e-2 5*1e1])
ylim([7*1e-3 5*1e0])
