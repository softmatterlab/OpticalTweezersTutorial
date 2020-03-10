function [fc_exp,D_exp,Efc_exp,ED_exp,f,XX,fw_mean,Pk,fcut,h]=psdfit_analytic(Xs,dt,Nb,pc)
% PSDFIT implementation of the PSD calibration method
%
%[fc_exp,D_exp,Efc_exp,ED_exp,h]=PSDFIT(Xs,dt,Nb)
% Efc and ED are useful if there is only one series of time.
% Nb: Number of windows, for example Nb=500
% pc: fraction of f_Niquist taken for the fitting. Usually 1/4 is a good value. 
[N,Nexp]=size(Xs);

Ts=N*dt;

f=1/dt*(1:(N/2))/N;

fwin=linspace(f(1),f(end),Nb);

fw_mean=(fwin(1:end-1)+fwin(2:end))/2;

for j=1:Nexp
    
    FT=fft(Xs(:,j));
    
    xx=FT.*conj(FT)*dt.^2/(N*dt);
    
    XX=xx(2:N/2+1);  
    
    %windowing
    
    Pk=[];
    
    nw=[];
    
    for k=1:length(fwin)-1
    
        ind=find(f>=fwin(k)&f<fwin(k+1));
        
        Pk(k)=mean(XX(ind));
        
        nw(k)=length(ind);
    end
    
    ew_inv=(1./Pk)./sqrt(nw);
    
    % Cropping the data
    fcut=pc*1/(2*dt);
    
    indf=find(fw_mean<=fcut);
    
    fwc=fw_mean(indf);
    
    Pkc=Pk(indf);
    
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
    
    fc_exp(j)=sqrt((S(1,2)*S(3,3)-S(2,2)*S(2,3))/(S(2,2)*S(1,3)-S(1,2)*S(2,3)));
    
    D_exp(j)=2*pi^2*(S(1,3)*S(3,3)-S(2,3)^2)/(S(2,2)*S(1,3)-S(1,2)*S(2,3));
    
    %k_exp(j)=gamma*2*pi*fc_exp(j);
    
    % Errors
    
    xmin=min(fwc)/fc_exp(j);
    
    xmax=max(fwc)/fc_exp(j);
    
    u=2*xmax/(1+xmax^2)-2*xmin/(1+xmin^2)+2*atan((xmax-xmin)/(1+xmin*xmax));
    
    v=4/(xmax-xmin)*atan((xmax-xmin)/(1+xmin*xmax))^2;
    
    sfc=sqrt(pi/(u-v));
    
    sD=sqrt(u/((1+pi/2)*(xmax-xmin)));
    
    Efc_exp(j)=sfc/sqrt(pi*fc_exp(j)*Ts)*fc_exp(j);
    
    ED_exp(j)=sqrt((1+pi/2)/(pi*fc_exp(j)*Ts))*sD*D_exp(j);
    
    %Ek_exp(j)=2*pi*gamma*Efc_exp(j);
    
end

%plot the last experiment

colors=[0,0.6,0.8];

h=figure(1);

clf

set(gcf,'Position',[150 300 800 600])

%Plot the PSD in loglog scale
axes('Position',[0.13 0.13 0.99-0.13 0.99-0.13])

loglog(f,XX*1e12,'.','MarkerSize',0.5,'Color',colors)

hold on

loglog(fw_mean,Pk*1e12,'.k','MarkerSize',6)

loglog(f,D_exp(end)/(2*pi^2)./(fc_exp(1)^2+f.^2)*1e12,'r','LineWidth',1)

loglog(fcut*ones(1,300),exp(linspace(log(0.8*min(XX)*1e12),log(1.1*max(XX)*1e12),300)),'.k','MarkerSize',2)

set(gca,'FontSize',16)

axis([min(f),max(f),0.9*min(XX)*1e12,1.1*max(XX)*1e12])

xlabel('$f_k(Hz)$','Interpreter','Latex')

ylabel('$|\hat{x}|^2/T_s, \, P^{(ex)}_k(\mu m^2/Hz)$','Interpreter','Latex')

hold off
