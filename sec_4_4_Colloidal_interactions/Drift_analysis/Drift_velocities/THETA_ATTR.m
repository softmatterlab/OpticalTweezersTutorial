function [ THETA ] = THETA_ATTR( x )
%THETA_ATTR: Scaling function for critical Casimir potential (attractive)
%   Input parameters:
%           x : l/csi
%
%   Output parameters:
%           THETA : scaling function

% Q6_coeff = [ -2.57611  -3.29886   2.27325   2.005*1e-1  1.0743*1e-2  -1.8317*1e-3  7.2396*1e-5];
% Log_coeff = [ 0.628319 -1.74764 ];
%
%
% if x>6
%
%     THETA = -9.425*x.*exp(-x);
%
% else
%
%     THETA = Q6_coeff*[ones(size(x)); x; x.^2; x.^3; x.^4; x.^5; x.^6]+ ( Log_coeff*[ones(size(x)); x]).*(x.*log(x));
%
% end


pQ6 = [ 7.2396*1e-5  -1.8317*1e-3  1.0743*1e-2  2.005*1e-1  2.27325  -3.29886  -2.57611 ];
pxLogx = [-1.74764 0.628319];


if x>6
    
    THETA = -9.425*x.*exp(-x);
    
else
    
    THETA = polyval(pQ6,x)+ polyval(pxLogx,x).*(x.*log(x));
    
end


end

