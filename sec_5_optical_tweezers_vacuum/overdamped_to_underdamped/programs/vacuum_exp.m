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


%analysis for the experimental results of vacuum trapping



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


%load data file with the psds, frequencies, pressures and temperatures

load('rawdata.mat')



pressure=pressures(2);
temperature=temperatures(2);
psd=psds(2, :);
frequency=frequencies(2,:);
x0=1;
x = lsqcurvefit(lorentz,x0,xdata,ydata)

