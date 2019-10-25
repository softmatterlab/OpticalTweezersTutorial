clear all 
close all

load('Data_positions_Fig9_1P2_S.mat')
kB=1.38e-23; % Boltzmann constant [m^2kg/s^2K]

r=1.03*1e-6;      % Particle radius [m]
v=0.00002414*10^(247.8/(-140+T));  % Water viscosity [Pa*s]
gamma_0=pi*6*r*v;
kth=1.28*1e-5;
xx = reshape(x,[size(x,1)*size(x,2),1]);
N=length(xx);
s_xn2=sum(xx.^2);

C_0=1/gamma_0+dt*s_xn2;
gamma_N=1/C_0;
sxn1xn=sum(xx(1:end-1).*xx(2:end));
K_0=kth/gamma_0;
K_N=1/dt*(-1+(dt*sxn1xn+(1+K_0*dt)/gamma_0)/(C_0));
a_0=1;
B_0=1e-6;
a_N=a_0+N/2;

B_N=B_0+1/2*sum((xx(2:end)-xx(1:end-1)*(1+K_N)).^2)/dt+1/(2*gamma_0)*(K_N-K_0).^2;
Dex=B_N/(2*(a_N-1));

kx=K_N*gamma_0;
