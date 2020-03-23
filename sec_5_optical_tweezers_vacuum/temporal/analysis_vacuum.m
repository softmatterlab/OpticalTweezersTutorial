%%code for vacuum trapping analyisis


%definitions
kB = 1.38064852e-23;
hbar = 1.0545718001391127e-34;
N_A = 6.022140857e+23;
c = 299792458.0;
eps0 =8.854187817620389e-12;

M_air = 28.97e-3; % kg/mol
T0 = 300; % Kelvins
viscosity = 18.27e-6; % Pascal seconds
mgas = M_air / N_A; % kg





% default values
radius = 136e-9/2; 
density = 1850;
refractive_index = 1.45+1j*6e-8;
nbb = 1.5+0.1*1i;
lambda_L = 1064e-9;
c_dp=2/5;
c_acc = 0.61;
gamma_sh = 7/5;


% Optics
wavelength = 1064e-9;
focal_length = 4e-3;
NA = 0.8;
aperture_radius_real = 7.5e-3/2;
width_inc_real = 1.06e-3 * 100/30;


P0 = 70e-3 ; % precisely the power that passes through the aperture
I_opt = P0*(2*pi/lambda_L)^2*NA^2 / (2*pi);% intensity at laser focus



eps_bb = nbb^2;

imag((eps_bb-1) /  (eps_bb+1)) , eps_bb;




pressure = logspace(-9,3, 100) * 100; % pressure in mbar from 1e-7 to  about atmosphere


T_int = T_internal(pressure, radius, refractive_index, I_opt, lambda_L, T0, nbb, c_acc, gamma_sh);
gamma_FB = 20*2*pi * ones(size(pressure)); % 100Hz damping from fb

[gamma_rad, S_rad] = damping_radiation(c_dp, I_opt, density, refractive_index, radius, lambda_L);

[gamma_gas, S_gas] = damping_gas(pressure);
gamma_tot = gamma_gas + gamma_rad + gamma_FB;
S_noise = S_gas + S_rad;
Gamma = pi*S_noise / particle_mass();
T_FB = Gamma/gamma_tot / kB;

gamma_gas_2_bath, S_gas_2_bath = damping_gas_2_bath(pressure, T_int);
gamma_tot_2_bath = gamma_gas_2_bath + gamma_rad + gamma_FB;
S_noise_2_bath = S_gas_2_bath + S_rad;
Gamma_2_bath = pi*S_noise_2_bath / particle_mass();
T_FB_2_bath = Gamma_2_bath/gamma_tot_2_bath / kB;

function mfp=mean_free_path(Pgas) 
mfp=sqrt(pi*kB*T0 / (2*mgas))*viscosity / Pgas;
end

function particle_mass=particle_mass(radius, density)

    volume(radius)*density;
end

function volume

volume=@(radius) 4*pi/3*radius^3;

end


function [gamma, damping]= damping_gas(Pgas,T, radius, density)
%gas damping and force sepctral density from low to high viscosiy regime when Tint = Tgas
%Pgas = gas pressure in Pascal
% returns gas damping in 2*pi *  Hertz (gamma)

m = particle_mass(radius, density);
mfp = mean_free_path(Pgas);
Kn = mfp / radius;
cK = 0.31*Kn / (0.785+1.152*Kn+Kn^2);
gamma = 2*np.pi * 3*radius*viscosity/m*0.619/(0.619+Kn)*(1+cK)
damping=kB * T * m*gamma / pi;

end

function [damping, fsd]=damping_gas_2_bath(Pgas, radius, density, Tgas, Tint ,c_acc)
% gas damping and force sepctral density in low viscosiy regime where Tint can be different from Tgas
% 
% Pgas = gas pressure in Pascal
% 
% returns gas damping in 2*pi *  Hertz (gamma)

Tem = c_acc * (Tint-Tgas)+Tgas;
m = particle_mass(radius, density);

gamma_im  = np.sqrt(8*mgas / (pi * kB * Tgas)) * Pgas / (density * radius);
gamma_em = pi/8 * np.sqrt(Tem / Tgas) * gamma_im;

damping = gamma_im + gamma_em;
fsd= m/pi * kB * (gamma_im * Tgas + gamma_em* Tem);
end

function polarizability=polarizability(refractive_index, radius , lambda_L)

eps_p  =refractive_index^2;
k_L = 2*pi/lambda_L;
V = volume(radius);
chi = 3 * (eps_p-1) / (eps_p+2);
a0 = eps0 * V * chi;
polarizability=a0 / (1-1i * k_L^3*a0 / (6*pi*eps0));
end

function [damping_radiation, noise]= damping_radiation(c_dp, I_opt, density, refractive_index , radius, lambda_L)
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


function power_absorbed=power_absorbed(refractive_index, radius, I_opt, lambda_L)
 
V = volume(radius)
eps_p =refractive_index^2
power_absorbed= 6*pi*I_opt/lambda_L* V * np.imag((eps_p-1)/(eps_p+2))
end


function power_convective=power_convective(Pgas,T_bulk , c_acc, radius,gamma_sh, Tgas)
    
%returns convective power and the first derivative wrt temperature
Tbulk=T0; 
c_acc=c_acc;
radius=radius;
Tgas=T0;
vrms = np.sqrt(3*kB*Tgas / mgas);
power_convective=c_acc*sqrt(2/(3*pi))*(pi*radius)^2*Pgas*vrms *(gamma_sh+1) / (gamma_sh-1)*(Tbulk / Tgas - 1);

end


function power_blackbody= power_blackbody(T, nbb, radius)

returns radiative power due to black body radiation and derivative wrt temperature
eps_bb = nbb^2;
zeta_R = 1.04; % rieman zeta (5)
V = volume(radius);
power_blackbody=72*zeta_R/pi^2 * V / (c^3*hbar^4)*np.imag((eps_bb-1) /  (eps_bb+1)) * (kB*T)^5;
end
    

function T_int=T_internal(Pgas, radius, refractive_index, I_opt, lambda_L,T0, nbb, c_acc, gamma_sh)

%Pgas: pressure in Pascal (default 101325 Pa= 1 atm ~ 1Bar)
 
    
    
    function f=func(T)
        Pconv = -power_convective(Pgas, T,c_acc, radius, gamma_sh, T0);
        Pabs = power_absorbed(refractive_index, radius, I_opt, lambda_L);
        Pbb_abs = power_blackbody(T0, nbb, radius);
        Pbb_em = -power_blackbody(T , nbb, radius);
       f=Pconv+Pabs+Pbb_abs+Pbb_em;
    end


    % check if input is array
%     if len(size([Pgas])) == 2: 
%         is_input_array = True
%     else:
%         is_input_array = False
%         Pgas = np.array([Pgas])
%         

    T_int = zeros(size(Pgas));
    
    for i=1:length(Pgas)
        P=i;
        pressure = P;
        
        sol = py.optimize.bisect(func, T0*0.9, 10*T0);
        T_int(i) = sol;
    end
end 



