function [fc_exp,D_exp,Efc_exp,ED_exp,f,XX,fw_mean,Pk,EPk,fcut]=psd_lfit(Xs,dt,Nb,pc)
% PSD_LFIT implementation of the PSD calibration using linear fitting
%
%[fc_exp,D_exp,Efc_exp,ED_exp,h]=PSD LFIT(Xs,dt,Nb)
% Efc and ED are useful if there is only one series of time.
% Nb: Number of windows, for example Nb=500
% pc: fraction of f_Niquist taken for the fitting. Usually 1/4 is a good value. 
[N,Nexp]=size(Xs);

kb=1.38064852e-23;
Ts=N*dt;

f=1/dt*(1:(N/2))/N;

fwin=linspace(f(1),f(end),Nb);

fw_mean=(fwin(1:end-1)+fwin(2:end))/2;

% Cropping the frequencies
fcut=pc*1/(2*dt);

indf=find(fw_mean<=fcut);

fwc=fw_mean(indf);

for j=1:Nexp
    
    x=Xs(:,j);
    
    FT=fft(x);
    
    xx=FT.*conj(FT)*dt.^2/(N*dt);
    
    XX=xx(2:N/2+1);  
    
    %windowing
    
    Pkw=[];
    
    nw=[];
    
    for k=1:length(fwin)-1
    
        ind=find(f>=fwin(k)&f<fwin(k+1));
        
        Pkw(k)=mean(XX(ind));
        
        nw(k)=length(ind);
    end
    
    Pkj(:,j)=Pkw';
    
end

Pk=mean(Pkj,2)';

EPk=Pk./sqrt(nw*Nexp);

Pkc=Pk(indf);

EPkc=EPk(indf);

%eliminate NaNs

ind=find(isnan(Pkc));

Pkc(ind)=[];

fwc(ind)=[];

% Fitting by analytical equations

for p=1:3
    
    for q=1:3
        
        S(p,q)=sum(fwc.^(2*(p-1)).*Pkc.^(q-1));
        
    end
    
end

%Parameters

fc_exp=sqrt((S(1,2)*S(3,3)-S(2,2)*S(2,3))/(S(2,2)*S(1,3)-S(1,2)*S(2,3)));

D_exp=2*pi^2*(S(1,3)*S(3,3)-S(2,3)^2)/(S(2,2)*S(1,3)-S(1,2)*S(2,3));

% Errors

xmin=min(fwc)/fc_exp;

xmax=max(fwc)/fc_exp;

u=2*xmax/(1+xmax^2)-2*xmin/(1+xmin^2)+2*atan((xmax-xmin)/(1+xmin*xmax));

v=4/(xmax-xmin)*atan((xmax-xmin)/(1+xmin*xmax))^2;

sfc=sqrt(pi/(u-v));

sD=sqrt(u/((1+pi/2)*(xmax-xmin)));

Efc_exp=sfc/sqrt(pi*fc_exp*Ts)*fc_exp;

ED_exp=sqrt((1+pi/2)/(pi*fc_exp*Ts))*sD*D_exp;






end
