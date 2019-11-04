function [k_psd, Ek_psd, mgamma_psd, Egamma_psd]=plotsub_msd(filename, positioninthefig1, title1, T, subs, maxlag)
load(filename);
disp(filename);
kb=1.38e-23;

col3=[0.00,0.45,0.74];

axes( 'Position',positioninthefig1);  % fa in modo di centrare il riquadro degli assi nella posizione voluta


% PSD analysis for Nexp

% Initialization of the workspace

addpath psd

nw=round(size(x(1:subs:end,:),1)/500); 
kb=1.38064852e-23;
gamma=6*pi*eta*a;
[fc_exp,D_exp,Efc_exp,ED_exp,f,XX,fw_mean,Pk,EPk,fcut]=psd_lfit(x(1:subs:end,:),dt*subs,nw,1/4);
mgamma_psd=kb*T./D_exp;

Egamma_psd=kb*T./D_exp^2*ED_exp;

k_psd=2*pi*gamma*fc_exp;

Ek_psd=2*pi*gamma*Efc_exp;

% estimation of k using the estimated gamma
%mk2_psd=2*pi*mgamma_psd.*mfc_psd;

%Ek2_psd=2*pi*mgamma_psd.*Efc_psd+2*pi*mfc_psd*Egamma_psd;


plot(f,XX*1e18,'.','MarkerSize',6,'Color',col3,'DisplayName', 'Experimental power spectral density')

hold on

plot(fw_mean,Pk*1e18,'Color', 'k'  ,'MarkerSize',10, 'MarkerEdgeColor','k','DisplayName','Mean of the experimental power spectral density', 'LineWidth', 3)

plot(f,D_exp/(2*pi^2)./(fc_exp^2+f.^2)*1e18,'--', 'LineWidth',3,'Color','r', 'DisplayName', 'Linear fit')

plot(fcut*ones(1,300),exp(linspace(log(0.8*min(XX)*1e18),log(1.1*max(XX)*1e21),300)),'--k','MarkerSize',2, 'HandleVisibility', 'off')
text(2e3,1,'$f_{cut}$','Interpreter','latex','FontSize',20)
%text(tau0_exp_lf*ntaus,0.1e-4,[num2str(ntaus),'$\tau_0$'],'Interpreter','latex','FontSize',20)
% set(gca,'FontSize',16, )

xlabel('$f_k(\rm Hz)$','Interpreter','Latex',  'FontSize',30)
ylim([1e-9 5e1])
xlim([1e-2 5e4])

  xticks([1e-2,1e-1,1e0,1e1,1e2,1e3,1e4]);
  yticks([1e-9,1e-6,1e-3,1e0,1e3]);
 
%  set(gca, 'YTick', [1e-9: 1e-1 : 1e1])
 
ylabel('$|\hat{x}|^2/T_s, \, P^{(ex)}_k(\rm nm^2/Hz)$','Interpreter','Latex', 'FontSize',30)
box on, 
% set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'xscale', 'log','yscale', 'log', 'XMinorTick', 'on','YMinorTick', 'on', 'TickLength',[0.02, 0.01], 'XTickLabel',{'0.01','','1','','100','','10000'});
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'xscale', 'log','yscale', 'log', 'XMinorTick', 'on','YMinorTick', 'on', 'TickLength',[0.02, 0.01]);
%msd_fit= 2*kb*T/k_msd*(1- exp(-tau/tau0));
%plot(tau, msd_fit*1e12, 'LineWidth',3,'Color',col1, 'DisplayName',  'Non linear fitting');
%hold on 
%ntaus=6;
%errorbar(tau,  mmsd*1e12, Emsd*1e12,'o','MarkerSize',7 ,'LineWidth', 1.5,'Color',colbar, 'DisplayName', 'Experimental mean square displacement');
%box on
%xticks((-0.5:0.1:0.5)*1e-7);
%xlim([0 0.009]);
%ylim([0, 6.5]*1e-4)
%set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',15);
%xlabel('$\tau(\rm s)$','Interpreter','Latex', 'FontSize',20)
%ylabel('$\rm MSD (\rm{\mu m^{2}})$','Interpreter','Latex', 'FontSize',20)
%plot([tau0*ntaus,tau0*ntaus],[0,6.5*1e-4],'--k', 'HandleVisibility','off')
%text(tau0*ntaus,0.05*max(mmsd)*1e12,[num2str(ntaus),'$\tau_0$'],'Interpreter','latex','FontSize',20)

hold off
% title(title1)
% legend

xwi = 400;    % width of the plot square
bx1 = 150;     % extra space at the left
bx2 = 30;     % extra space at the right

Xpix = 3*xwi+3*bx1+3*bx2;  % total

ywi = 300;    % length riquadro con funzione
by1 = 110;     % extra space below
by2 = 50;     % extra space up

Ypix = 1*by1+1*ywi+1*by2;  % larghezza figura in pixel
 axes('Position',[(0) 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
hold on

xt = [bx1-70,bx1+xwi+bx1-45,2*bx1+2*xwi+bx1-5];
yt = [ by1+ywi+20,by1+ywi+20,by1+ywi+20];
str = {'\bf a','\bf b','\bf c'};
text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])

end