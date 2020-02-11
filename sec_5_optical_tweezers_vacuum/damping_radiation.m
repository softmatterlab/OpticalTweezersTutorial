function [damping_rad, noise]= damping_radiation(c_dp, I_opt, density, refractive_index , radius, lambda_L)
% radiation damping and noise due to photon shot noise

alpha = polarizability(refractive_index, radius, lambda_L);
m = particle_mass(radius, density);
k_L = 2*pi/lambda_L;
omega_L = k_L*c;
s_scatt = abs(alpha)^2 * k_L^4 / (6*pi*eps0^2);
Pscatt = s_scatt * I_opt;

print(sqrt(s_scatt));
damping_rad =c_dp*Pscatt / (m*c^2);
noise=c_dp*hbar*omega_L*Pscatt / (2*pi*c^2);
end