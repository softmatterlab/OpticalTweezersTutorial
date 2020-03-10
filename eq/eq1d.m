function [k_eq, sigma2_k_eq]=eq1d(x,T)
%function [mk, Ek]=eq1d(Vx,T,deltax)
% deltax: error in the position detection
% EQ1D   1D implementation of the EQUIPARTITION METHOD

kb=1.38064852e-23;

[N,Nexp]=size(x);
deltax=0;
for j=1:Nexp
    
    xx=x(:,j);
    
    sigma2=var(xx);    
    
    k=kb*T./(var(xx,[],1)-deltax^2);
    
    sigma2_kexp2(j)=4/(Nexp*N)*k^2/sigma2*deltax^2;
    
    kexp(j)=k;
    
end
k_eq=mean(kexp);

sigma2_k_eq=sqrt(var(kexp)+mean(sigma2_kexp2));

%
%disp('...')

%disp('Equipartition analysis')

%disp(['k_eq: ' num2str(k_eq) '+-' num2str(sigma2_k_eq)]);