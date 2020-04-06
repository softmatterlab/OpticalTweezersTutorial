function [gamma, damping]= damping_gas(Pgas,T, radius, density, viscosity)
%gas damping and force sepctral density from low to high viscosiy regime when Tint = Tgas
%Pgas = gas pressure in Pascal
% returns gas damping in 2*pi *  Hertz (gamma)
kB= 1.38064852*1e-23;
m = particle_mass(radius, density);
mfp = mean_free_path(Pgas, m,viscosity,T);
Kn = mfp / radius;
cK = 0.31*Kn / (0.785+1.152*Kn+Kn^2);
gamma = 2*pi * 3*radius*viscosity/m*0.619/(0.619+Kn)*(1+cK);
damping=kB * T * m*gamma / pi;

end