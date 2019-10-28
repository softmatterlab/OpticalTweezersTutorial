
clear all 
close all
%load data files
load('Data_positions_Fig9_1P4_S.mat')

xx = reshape(x,[size(x,1)*size(x,2),1]);
N=length(xx);
dnn=1e2*(1:1e3:50000);
for nj=1:length(dnn)
xn =xx(1:nj);
N=length(xn);
%some useful calculations from the positions
sxn1xn=sum(xn(1:end-1).*xn(2:end));
s_xn2=sum(xn.^2);

%physical constants
kB=1.38e-23; % Boltzmann constant [m^2kg/s^2K]

%definition of parameter

eta=0.00002414*10^(247.8/(-140+T));  % Water viscosity [Pa*s]
k_th=14*1e-6;
gamma_th=6*pi*eta*a;
D_th=kB*T/gamma_th;
sigma2_D_th=D_th/10;
sigma2_F_k_ex_th_gamma_th=k_th/(10*gamma_th);

%hyperparameters
psi_0=sigma2_F_k_ex_th_gamma_th/2*D_th;
beta_0=2*D_th*(D_th^2+sigma2_D_th)/sigma2_D_th;
K_0=k_th/gamma_th;

alpha_0=2+D_th^2/sigma2_D_th;


%parameters
denGN_0=1/psi_0+dt*s_xn2;

psi_N=1/denGN_0;

alpha_N=alpha_0+N/2;

K_N=1/dt*(-1+(dt*sxn1xn+(1+K_0*dt)/psi_0)/denGN_0);

beta_N=beta_0+1/2*sum((xn(2:end)-xn(1:end-1)*(1+K_N*dt)).^2)/dt+1/(2*psi_0)*(K_N-K_0).^2;


Dex=beta_N/(2*(alpha_N-1));
gamma(nj)=kB*T/Dex;

sigma2_D_ex=beta_N^2/(4*(alpha_N-1)^2*(alpha_N-2));
k_ex(nj)=K_N*gamma(nj);
sigma2_f_k_x_ex_gamma=psi_N*beta_N/(alpha_N-1);

end

figure(1)
plot(dnn, k_ex)
figure(2)
plot(dnn, gamma)
