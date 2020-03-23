function lorentz_overdamped=lorentz_overdamped(pressure,  density, temperature, g0, viscosity,coef, f) %lorentz_overdampded( pressure,  density, temperature, g0, viscosity,c, f0,  radius, f)
kB= 1.38064852*1e-23;
mass=particle_mass(coef(3), density);
if isempty(g0)
    [g0, ~] = damping_gas(pressure*100,temperature, coef(3), density, viscosity);
    g0 =g0/(2*pi);
end
fc = coef(2).^2./g0;
lorentz_overdamped=1/(2*pi.^3) * coef(1).^2 * kB * temperature ./ g0 ./ mass ./ (fc.^2+f.^2);
end
    