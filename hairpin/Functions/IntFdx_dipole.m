function [ Energy ] = IntFdx_dipole(d,force,KbT)
% Gives the integral from 0 force to force using the extension provided by dipole.
% 
    f_min= 0.000001;
    f_max = force;
	
    x_min= dipole(d,f_min,KbT)
    x_max= dipole(d,f_max,KbT)

    %Integrating in f:
    Energyf= d*(log(sinh(f_max*d/KbT))-log(f_max*d/KbT))-d*(log(sinh(f_min*d/KbT))-log(f_min*d/KbT))
    %Integrating in x instead of in f:
    Energy = (f_max-f_min)*(x_max-x_min)-Energyf;
end
