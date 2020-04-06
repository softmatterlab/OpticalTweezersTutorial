function [ Work ] = IntFdx_WLC(fmax,P,L,KbT)
% Gives the integral taking X as the integration variable and f_WLC(x) as
% the function to integrate! The result is pN*nm
% Courtesy of Aurelien Severino

x_min = 0;
x_max = L * WLC_inv(fmax,P, KbT);
x_int = linspace(x_min,x_max,1000);

Work = trapz(x_int, WLC(x_int./L, P, KbT, 0));

 
 
end

