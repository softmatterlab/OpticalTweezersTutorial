function [k, sigma_k, gamma, sigma_gamma, D, sigma_D]= bayesian(xn, dt,T, a, nsubs)
dt=dt*nsubs;
N=floor(length(xn)/nsubs);
s2=sum(xn(1:nsubs:end-nsubs).*xn(1+nsubs:nsubs:end));
s1=sum(xn(1:nsubs:end-nsubs).^2);
s3=sum(xn(1+nsubs:nsubs:end).^2);
%physical constants
kB=1.38064852e-23; % Boltzmann constant [m^2kg/s^2K]

%definition of parameter

eta=0.00002414*10^(247.8/(-140+T));  % Water viscosity [Pa*s]
k_th=14*1e-6;
gamma_th=6*pi*eta*a;
D_th=kB*T/gamma_th;
sigma2_D_th=(D_th/10).^2;
sigma2_F_k_ex_th_gamma_th=(k_th/(10*gamma_th)).^2;

%hyperparameters
psi_0=sigma2_F_k_ex_th_gamma_th/(2*D_th);
beta_0=2*D_th*(D_th^2+sigma2_D_th)/sigma2_D_th;
K_0=k_th/gamma_th;

alpha_0=2+D_th^2/sigma2_D_th;


%parameters
denGN_0=1/psi_0+dt*s1;%

psi_N=1/denGN_0;%

alpha_N=alpha_0+N/2;%

K_N=-1/dt*(-1+(dt*s2+(1-K_0*dt)/psi_0)/denGN_0);%
beta_N=beta_0+1/(2*psi_0)*(K_N-K_0)^2+(s3-2*s2*(1-K_N*dt)+s1*(1-K_N*dt)^2)/(2*dt);

D=beta_N/(2*(alpha_N-1));
gamma=kB*T/D;

sigma_D=sqrt(beta_N^2/(4*(alpha_N-1)^2*(alpha_N-2)));
sigma_gamma=kB*T/D^2*sigma_D;
k=K_N*gamma;
sigma_f_k_gamma=sqrt(psi_N*beta_N/(alpha_N-1));
sigma_k=gamma*sigma_f_k_gamma+K_N*sigma_gamma;
end