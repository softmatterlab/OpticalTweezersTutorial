function [z] = WLC_invFast2(f,P,KbT)
% [z] = WLC_inv2(f,P,KbT)
% P is the persistence length in [nm], KbT is in [pN*nm], the force is
% given in [pN]. The force must be positive.
% IMPORTANT: the function returns the value "z" of the normalized extension
% z=x_stretch/contour_length. Hence z is always in the interval [0,1]
% NOTE 1: this function works with f as a vector or scalar argument.
% Optionnaly, KbT can also be provided as vector of a size identical to f.
% Function courtesy of Aurélien Severino

if length(KbT)~=1 && length(KbT)~=length(f)
   error('KbT must be either of size 1 or of size n=length(f)')
end

if any(f<0)
    disp('Negative forces are not allowed')
end

fnorm = ((4*P)./KbT).*f;
z = zeros(size(f));

for j = 1:length(f)
    
    a2 = (1/4)*(-9-fnorm(j));
    a1 = (3/2)+(1/2)*fnorm(j);
    a0 = -fnorm(j)/4;
    
    R = (9*a1*a2-27*a0-2*a2^3)/54;
    Q = (3*a1-a2^2)/9;
    
    D = Q^3+R^2; % the generalised discriminant for cubic equations
    
    if D > 0
        % In this case, there is only one real root, given by "out" below
        S   = nthroot(R+sqrt(D),3);
        T   = nthroot(R-sqrt(D),3);
        out = (-1/3)*a2+S+T;
    elseif D < 0
        % In this case there 3 real distinct solutions, given by out1,
        % out2, out3 below. The one that interests us is that in the
        % inerval [0,1]. It is seen ("empirically") that is always the
        % second one in the list below [there is perhaps more to search here]
        
        theta = acos(R/sqrt(-Q^3));
        %out1 = 2*sqrt(-Q)*cos(theta/3)-(1/3)*a2;
        out2 = 2*sqrt(-Q)*cos((theta+2*pi)/3)-(1/3)*a2;
        %out3 = 2*sqrt(-Q)*cos((theta+4*pi)/3)-(1/3)*a2;
        
        % We implement the following check just to be sure out2 is the good root 
        % (in case this "empirical" truth turns out to stop working) 
        if out2 < 0 || out2 > 1 
            error(['The default root doesn"t seem the be good one - you ',...
                'may want to check if the others lie in the interval [0,1]'])
        else
            out = out2;
        end
    else
        % In theory we always go from D>0 to D<0 by passing to a D=0
        % boundary, where we have two real roots (and where the formulas
        % above change again slightly). In practice, however, due to round-off errors,
        % it seems we never hit this boundary but always pass "through" it 
        % This D=0 scenario could still be implemented if needed, though.
        error('Lalala')
    end
    z(j) = out;
end

    
end
