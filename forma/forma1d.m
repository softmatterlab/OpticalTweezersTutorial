function [fcv, Dv, Efc, ED] = forma1d(Vx, dt)
% FORMA1D   1D implementation of FORMA
%
% [fc, D] = FORMA1D(x, dt)
% Efc and ED are useful if there is only one series of time.

[N,Nexp]=size(Vx);

for j=1:Nexp
    
    x=Vx(:,j);
    
    dx = diff(x)/dt;
    
    x = x(1:end-1);
    
    fc= -1/(2*pi)*sum(x.*dx)/sum(x.^2);
    
    eps=dx+2*pi*fc*x;
    
    D=dt/2*mean((eps).^2);
    
    ED(j)=D*sqrt(2/N);
    
    Efc(j)=1/(2*pi)*sqrt(2*D/(dt*sum(x.^2)));
    
    fcv(j)=fc;
    
    Dv(j)=D;
    

    
end


