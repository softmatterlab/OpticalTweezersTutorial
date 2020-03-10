%This program computes the real and the imaginary part of the complex
%modulus of a viscoelastic material, measured by Active Microrheology.
%The particle radius and trap stiffness must be specified as an input
%A file containing the coordinates of the trapped article, the sinusoidal 
%perturbative force at a single frequency and the corresponding time, 
%is needed, and must be specified in the parameters 
%filename. 



clear all
close all
R = 1.0e-6; %particle radius
k = 1e-6; %trap stiffness (previously measure ny equipartition in thermal equilibrium)

%Write driving frequencies at which sinusoidal motion of the optial trap
%was applied, they must be specified in the file name
FF = [0.05 0.1 0.2 0.4 0.8 1.0 1.6 3.2 6.4 12.8 25.6];

Filnam ='ActiveMicrorheologyData/5mMdrive_k1e-06_freq_'

for i = 1:length(FF)

ff = FF(i);
    
filename = [Filnam  num2str(ff) '.mat']; 

load(filename)

%read data
fs=data.fs; %sampling frequency
x=data.x; %particle position
fd=data.fd; %perturbative sinusoidal force
t=data.t; %time

nfft = 2^17; %number of point to compute cross power spectral density

[pxx,f] = cpsd(x,x,[],[],nfft,fs);
[pxfd,f] = cpsd(fd,x,[],[],nfft,fs);

invchi = pxfd./pxx; %1/Fourier transform of response function of x to fd

reinvchi = real(invchi); %real part of 1/Fourier transform of response function of x to fd
iminvchi = imag(invchi); %imaginary part of 1/Fourier transform of response function of x to fd

figure(1)
loglog(f,pxx)
hold on
xlabel('f [Hz]')
ylabel('PSD [m^2/Hz]')
title('Power spectral density of x for particle subject to sinusoidal force at single frequency')

[mm,I]=(min(abs(ff - f))); 
freqd = f(I);

frqd(i) = freqd;
G1(i) = (reinvchi(I) - k)/6/pi/R %storage modulus
G2(i) = iminvchi(I)/6/pi/R %loss modulus

end

figure(2)
loglog(frqd,G1,'r o')
hold on
loglog(frqd,G2,'s b')
xlabel('f [Hz]')
ylabel('G(f),G''(f) [Pa]')
title('Storage (circles) and loss (squares) modulus of material obtained bt Active Microrheology')


