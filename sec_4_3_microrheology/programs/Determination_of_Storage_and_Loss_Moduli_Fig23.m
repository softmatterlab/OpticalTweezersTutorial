%% This program computes the real and the imaginary part of the complex
%modulus of a viscoelastic material, measured by Active Microrheology.
%The particle radius and trap stiffness must be specified as an input
%A file containing the coordinates of the trapped article, the sinusoidal 
%perturbative force at a single frequency and the corresponding time, 
%is needed, and must be specified in the parameters 
%filename. 



clear all
close all

col3=[0.00,0.45,0.74];
 
bx1 = 140;    % extra space at the left
xwi = 540;     % width of the plot square
bx2 = 20;      % extra space at the right

  
Xpix =1400;   %  width of figure in pixel
by1 = 95;      % extra space below
ywi = 450;    % length of the plot square
by2 = 50;     % extra space up


Ypix = by1+ywi+by2;  % width of the figure in pixel



%number of bins of the histogram, if not set default is 50
figure('Position',[10 20 Xpix Ypix]);
R = 1.0e-6; %particle radius
k = 1e-6; %trap stiffness (previously measure ny equipartition in thermal equilibrium)

%Write driving frequencies at which sinusoidal motion of the optial trap
%was applied, they must be specified in the file name
FF = [0.05 0.1 0.2 0.4 0.8 1.0 1.6 3.2 6.4 12.8 25.6];

Filnam ='../data/ActiveMicrorheologyData/5mMdrive_k1e-06_freq_'
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
if i==6
 xsel=x;
 tsel=t;
end
nfft = 2^17; %number of point to compute cross power spectral density

[pxx,f] = cpsd(x,x,[],[],nfft,fs);
[pxfd,f] = cpsd(fd,x,[],[],nfft,fs);

invchi = pxfd./pxx; %1/Fourier transform of response function of x to fd

reinvchi = real(invchi); %real part of 1/Fourier transform of response function of x to fd
iminvchi = imag(invchi); %imaginary part of 1/Fourier transform of response function of x to fd


loglog(f,pxx*1e18)
hold on

[mm,I]=(min(abs(ff - f))); 
freqd = f(I);

frqd(i) = freqd;
G1(i) = (reinvchi(I) - k)/6/pi/R; %storage modulus
G2(i) = iminvchi(I)/6/pi/R; %loss modulus

end

xlabel('$f \rm (Hz)$', 'Interpreter','Latex', 'FontSize',30)
ylabel('PSD (nm$^2$ Hz$^{-1}$)', 'Interpreter','Latex', 'FontSize',30)

xlim([1e-2 1e2])
 ylim([1e-2 1e7])

set(gca,'TickLabelInterpreter','Latex', 'linewidth',1.5,'FontSize',25, 'TickLength',[0.02, 0.01]);
hold off

positionintefig2=[bx1+100 0 220 0]/Xpix + [0 170 0 90]/Ypix;
axes( 'Position',positionintefig2);

plot(tsel, xsel*1e9)
hold on
plot(tsel,0.2*1e3*cos(2*pi*tsel),'LineWidth', 2,'Color','k')
xlim([0 4])
ylim([-220 220])

xlabel ('$t \rm (s)$', 'Interpreter','Latex', 'FontSize',30)
ylabel ('$x \rm (nm)$', 'Interpreter','Latex', 'FontSize',30)
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',18,'TickLength',[0.02, 0.01]);
%%
positionintefig3=[2*bx1+xwi+bx2 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix;
axes( 'Position',positionintefig3);
loglog(frqd,G1,'o', 'MarkerSize',10, 'Color','r','MarkerFaceColor','r')
hold on
loglog(frqd,G2,'s','MarkerSize',10, 'Color',col3,'MarkerFaceColor',col3)
xlabel('$f \rm (Hz)$', 'Interpreter','Latex', 'FontSize',30)
ylabel('$G''(f),G''''(f) \rm (Pa)$', 'Interpreter','Latex', 'FontSize',30)
%title('Storage (circles) and loss (squares) modulus of material obtained bt Active Microrheology')

set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01]);
%%


%we add the data of the passive case
NFFT = 2^17;
Filepath = '../data/PassiveMicrorheologyData/';
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




%Determine trapstifness by equipartition theorem
k=kB*T0/var(x) 

%Compute PSD
[f,spx] = spectave(detrend(x),fs,NFFT);



%Smooth PSD using polynomial fit of degree D
D =11;
logsx=log(spx(3:end));
logf = log(f(3:end));

coeffx=polyfit(logf,logsx,D);

spx1=exp(polyval(coeffx,(logf),D));
spx(3:end)=spx1;



%Compute imaginary (chi2) and real (chi1) parts of linear response function of x
%to perturbative force
chi2=pi*f.*spx/2/kB/T0; %Fluctuation-dissi�tion Theorem
chi11=KramersKronig(2*pi*f,chi2');
chi1=chi11';

%Compute storage (G1) and loss (G2) modulus from linear response function
G1 = chi1./(chi1.^2+chi2.^2)/6/pi/r - k/6/pi/r;
G2 = chi2./(chi1.^2+chi2.^2)/6/pi/r;



loglog(f(4:2500),G1(4:2500),'r --','LineWidth', 3 )

loglog(f(4:2500),G2(4:2500),'--', 'Color',col3,'LineWidth', 3)

xlim([3e-2 4e1])
ylim([3e-3 2e0])
xticks([1e-1,1e0,1e1]);
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01]);


axes('Position',[(0) 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix);  
hold on

xt = [bx1-100, bx1+xwi+bx2+20];
yt = [ by1+ywi+20, by1+ywi+20];
str = {'\bf a','\bf b'};
text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])


%saveas(gcf,'Fig23.eps','epsc')
%saveas(gcf,'Fig23.fig')