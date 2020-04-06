function lorentz=lorentz(c, f0,  radius, pressure,  density, temperature, g0, viscosity, f)
kB= 1.38064852*1e-23;
mass=particle_mass(radius, density);
if isempty(g0)
    [g0, ~] = damping_gas(pressure*100,temperature, radius, density, viscosity);
    g0 =g0/(2*pi);
end
lorentz=1./(2*pi^3) * c^2 * kB * temperature * g0 ./ mass ./ ((f0.^2-f.^2).^2 + g0.^2*f.^2);
end
    
