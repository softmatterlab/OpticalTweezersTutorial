function [force] = WLC(x_norm,P,KbT, bool)
% [force] = WLC(x_norm,P,KbT, bool)
% This functions compute the force in the WLC model using:
%   - bool = 0 the formula of Marko and Siggia (1995) 
%   - bool = 1 the Bouchiat interpolation improvement
%  Other input:  
%    - x_norm = extension / contour length
%    - P: persistence length
% 
% Note: this function works with x_norm as a vector or scalar argument, the rest must
% be scalar

if bool == 1
    force = (KbT./P).*(1./(4.*(1-x_norm).^2) -1/4 + x_norm ...
            - 0.5164228.*x_norm.^2  - 2.737418.*x_norm.^3 ...
            + 16.07497.*x_norm.^4 - 38.87607.*x_norm.^5 ...
            + 39.49944.*x_norm.^6 - 14.17718.*x_norm.^7);
elseif bool == 0
    force = (KbT./P) .* ( 1./(4.*(1-x_norm).^2) - 1/4 + x_norm );
else
    error('error in WLC.m')
end