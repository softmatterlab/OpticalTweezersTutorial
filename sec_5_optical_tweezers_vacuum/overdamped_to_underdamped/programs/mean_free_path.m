function mfp=mean_free_path(Pgas, mgas, viscosity, T0) 
kB = 1.38064852e-23;
mfp=sqrt(pi*kB*T0 / (2*mgas))*viscosity/ Pgas;
end