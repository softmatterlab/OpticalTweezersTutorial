function [extension] = dipole(d,force,KbT)
    % [extension] = dipole(d,force,KbT)
    x=(d.*force)./KbT;
    extension=d.*(coth(x)-1./x);
end
