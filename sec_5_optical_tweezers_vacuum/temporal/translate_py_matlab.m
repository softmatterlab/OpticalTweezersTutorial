function mfp=mean_free_path(Pgas)
kB = 1.38064852e-23;
hbar = 1.0545718001391127e-34;
N_A = 6.022140857e+23;
c = 299792458.0;
eps0 =8.854187817620389e-12;
   
mfp= sqrt(pi*kB*T0 / (2*mgas))*viscosity / Pgas;
end