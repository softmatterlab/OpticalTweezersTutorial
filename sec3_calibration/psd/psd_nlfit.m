function [fc_exp,D_exp,Efc_exp,ED_exp,f,XX,fw_mean,Pk,EPk,fcut,h]=psd_nlfit(Xs,dt,Nb,pc)
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

% Cropping the frquecies
fcut=pc*1/(2*dt);

indf=find(fw_mean<=fcut);

fwc=fw_mean(indf);

for j=1:Nexp
    
    FT=fft(Xs(:,j));
    
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

Pk=mean(Pkj,2);

EPk=Pk./sqrt(nw'*Nexp);

Pkc=Pk(indf);

EPkc=EPk(indf);

wc=1./EPkc.^2;

%eliminate NaNs

ind=find(isnan(Pkc));

Pkc(ind)=[];

fwc(ind)=[];

wc(ind)=[];

% First estimations to use in the nl fitting as starting points

a0=Pkc(end)*fwc(end)^2;

b0=a0/Pkc(1);

% Using nl fitting

ft=fittype('a/(x^2+b)');

maxpk=max(Pkc);

c=fit(fwc',Pkc/maxpk,ft,'Weights',wc'*maxpk^2,'StartPoint',[a0/maxpk,b0],'Display','Iter');

%Parameters

fc_exp=sqrt(c.b);

D_exp=2*pi^2*c.a*maxpk;

% Errors

cint=confint(c,0.95);

Efc_exp=1/2*1/(c.b)^(1/2)*(cint(2,2)-cint(1,2))/2;

ED_exp=2*pi^2*(cint(2,1)-cint(1,1))/2*maxpk;
    
%plot the last experiment

colors=[0,0.6,0.8];

h=figure(1);

clf

col1=[0.9,0.4,0.2];

set(gcf,'Position',[150 300 800 600])

%Plot the PSD in loglog scale
axes('Position',[0.13 0.13 0.99-0.13 0.99-0.13])

loglog(f,XX*1e12,'.','MarkerSize',0.5,'Color',colors)

hold on

loglog(fw_mean,Pkw*1e12,'ok','MarkerSize',5)

loglog(f,D_exp/(2*pi^2)./(fc_exp^2+f.^2)*1e12,'LineWidth',1,'Color',col1)

loglog(fcut*ones(1,300),exp(linspace(log(0.8*min(XX)*1e12),log(1.1*max(XX)*1e12),300)),'.k','MarkerSize',2)

set(gca,'FontSize',16)

axis([min(f),max(f),0.9*min(XX)*1e12,1.1*max(XX)*1e12])

xlabel('$f_k(Hz)$','Interpreter','Latex')

ylabel('$|\hat{x}|^2/T_s, \, P^{(ex)}_k(\mu m^2/Hz)$','Interpreter','Latex')

hold off
