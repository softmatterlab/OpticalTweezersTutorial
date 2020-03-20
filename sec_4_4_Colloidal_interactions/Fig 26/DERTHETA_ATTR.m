function [ DerTHETA ] = DERTHETA_ATTR( x )
%DERTHETA_ATTR: Scaling function for critical Casimir potential
%   Input parameters:
%           x : l/csi
%
%   Output parameters:
%           DerTHETA : scaling function
%
% Q6_coeff = [ -2.57611  -3.29886   2.27325   2.005*1e-1  1.0743*1e-2  -1.8317*1e-3  7.2396*1e-5];
% DQ6_coeff = [Q6_coeff(2:7).*[1 2 3 4 5 6] 0];
% Log_coeff = [ 0.628319 -1.74764 ];
%
% if x>6
%
%     DerTHETA = -9.425*(ones(size(x))-x).*exp(-x);
%
% else
%
%     DerTHETA = DQ6_coeff*[ones(size(x)); x; x.^2; x.^3; x.^4; x.^5; x.^6]+...
%         ( Log_coeff*[ones(size(x)); 2*x]).*log(x) + Log_coeff*[ones(size(x)); x];
%
% end

pQ6 = [ 7.2396*1e-5  -1.8317*1e-3  1.0743*1e-2  2.005*1e-1  2.27325  -3.29886  -2.57611 ];
pxLogx = [-1.74764 0.628319];

DpQ6 = pQ6(1:6).*(6:-1:1);
DpxLogx = pxLogx.*[2 1];

if x>6
    DerTHETA = -9.425*polyval([-1 1],x).*exp(-x);
else
    DerTHETA = polyval(DpQ6,x)+...
        polyval(DpxLogx,x).*log(x) + polyval(pxLogx,x);
end


end

