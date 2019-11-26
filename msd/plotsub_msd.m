function [k_msd,sigma_k_msd ,  gamma_msd, sigma_gamma_msd , D_msd , sigma_D_msd,tau0, sigma_tau0]=plotsub_msd(filename, positioninthefig1, title1, T, subs, maxlag,partau0,ytau,aa)
load(filename);
disp(filename);
kB=1.38064852e-23;
%blue color
col1=[73/255,4/255,10/255];
%yellow
col2=[241/255,185/255,14/255];
%gray color for experimental data
colbar=[7/255, 79/255, 129/255];


[k_msd,sigma_k_msd, tau0, sigma_tau0, D_msd, sigma_D_msd, tau, mmsd, Emsd, indc, gamma_msd, sigma_gamma_msd]=msd_nfilt(x(1:subs:size(x,1),:),T,dt*subs,maxlag);


% % plot
% figure(1)
% 
% clf
% 
% set(gcf,'Position',[150 300 800 600])
% 
% axes('OuterPosition',[0 0 1 1])
% 
% errorbar(tau(1:1:end),mmsd(1:1:end)*1e12,Emsd(1:1:end)*1e12,'ob','LineWidth',1)
% 
% hold on
% 
% plot(tau_cut(1:1:end),a0*(1-exp(-tau_cut(1:1:end)/tau0))*1e12,'b')
% 
% plot([ntaus*tau0,ntaus*tau0],[0.8*min(mmsd),1.3*max(mmsd)]*1e12,'--k')
% 
% axis([tau(1),tau(end),0.8*min(mmsd)*1e12,1.2*max(mmsd)*1e12])
% 
% text(ntaus*tau0,0.05*max(mmsd)*1e12,[num2str(ntaus),'$\tau_0$'],'Interpreter','latex','FontSize',16)
% 
% set(gca,'FontSize',16)
% 
% xlabel('$\tau$','Interpreter','latex')
% 
% ylabel('$MSD_x(\mu \textrm{m}^2)$','Interpreter','latex')
% 
% %   
% disp('...')
% 
% disp('MSD by non-linear fitting')
% 
% disp(['tau0_msd: ' num2str(tau0) '+-' num2str(Etau)])
% 
% disp(['k_msd: ' num2str(k_msd*1e6) '+-' num2str(Ek_msd*1e6)])
% 
% disp(['D_msd: ' num2str(D_msd) '+-' num2str(ED_msd)])
% 
% disp(['gamma_msd: ' num2str(gamma_msd) '+-' num2str(Egamma_msd)])
% 
col3=[0.00,0.45,0.74];
axes( 'Position',positioninthefig1);  % fa in modo di centrare il riquadro degli assi nella posizione voluta

msd_fit= 2*kB*T/k_msd*(1- exp(-tau/tau0));


ntaus=6;
e=errorbar(tau(1:2:end),  mmsd(1:2:end)*1e18, Emsd(1:2:end)*1e18,'.','MarkerSize',20,'LineWidth', 1.5, 'Color', colbar, 'DisplayName', 'Experimental MSD');
e.Color = col3;
hold on 
plot(tau, msd_fit*1e18, '--','LineWidth',3,'Color','r', 'DisplayName',  'Non -linear fitting');
box on
%xticks((-0.5:0.1:0.5)*1e-7);
xlim([0 0.009]);
ylim([0, 6.5]*1e2)



if aa==5
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5, 'FontSize',25);
ylabel('$\rm MSD (\rm{nm^{2}})$','Interpreter','Latex', 'FontSize',30)
else 
   set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5, 'FontSize',25,'Yticklabel',[]);
    end
 



xlabel('$\tau(\rm s)$','Interpreter','Latex', 'FontSize',30)





plot([tau0*ntaus,tau0*ntaus],[0,6.5*1e2],'--k', 'HandleVisibility','off')
% text(tau0*ntaus*partau0,0.05*max(mmsd)*1e12,[num2str(ntaus),'$\tau_0$'],'Interpreter','latex','FontSize',30)
text(tau0*ntaus*partau0,1.2*ytau*1E2,[num2str(ntaus),'$\tau_{{\rm ot},x}$'],'Interpreter','latex','FontSize',30)



 if aa==5
     
     LL= legend ({'Experimental MSD','Non-linear fitting'},'Box','off','Position',[0.18 0.25 0.1 0.2])

LL.FontSize = 18
 end
hold off
% title(title1)
% legend
xwi = 400;    % width of the plot square
bx1 = 120;     % extra space at the left
bx2 = 20;     % extra space at the right

Xpix = 3*xwi+bx1+3*bx2;  % total

ywi = 300;    % length riquadro con funzione
by1 = 110;     % extra space below
by2 = 70;     % extra space up

Ypix = 1*by1+1*ywi+1*by2;  % larghezza figura in pixel
 axes('Position',[(0) 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
hold on
% 
% xt = [bx1-100,bx1+xwi+bx1-75,2*bx1+2*xwi+bx1-40];
% yt = [ by1+ywi+20,by1+ywi+20,by1+ywi+20];


xt = [bx1-100,bx1+xwi+bx2,bx1+2*(xwi+bx2)];
yt = [ by1+ywi+30,by1+ywi+30,by1+ywi+30];
str = {'\bf a','\bf b','\bf c'};
text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])


end