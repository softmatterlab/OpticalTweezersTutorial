function [k_acf_lf,sigma_k_acf_lf,D_acf_lf, sigma_D_acf_lf,gamma_acf_lf, sigma_gamma_acf_lf , k_acf_nl, sigma_k_acf_nl, D_acf_nl, sigma_D_acf_nl,gamma_acf_nl, sigma_gamma_acf_nl, tau0_exp_lf, tau0_exp_nl,Xpix, Ypix]=plotsub_pot(filename, positioninthefig1, title1, T, subsample, partau0, end_plot,aa)
load(filename);
disp(filename);

 
%blue color

    
col3=[0.00,0.45,0.74];


colbar=[7/255, 79/255, 129/255];

[k_acf_lf, sigma_k_acf_lf, D_acf_lf, sigma_D_acf_lf,gamma_acf_lf, sigma_gamma_acf_lf,tau_acf_lf, mc, Ec,indc, tau0_exp_lf, c0_exp_lf]=acf_lfit(x(1:subsample:size(x,1),:),T,dt*subsample);

[k_acf_nl, sigma_k_acf_nl, D_acf_nl, sigma_D_acf_nl,gamma_acf_nl, sigma_gamma_acf_nl, tau_nl, mc, Ec, indc, tau0_exp_nl, sigma_tau0_nl,c0_exp_nl]=acf_nlfit(x(1:subsample:size(x,1),:),T,dt*subsample);
 
axes( 'Position',positioninthefig1);  % fa in modo di centrare il riquadro degli assi nella posizione voluta





e=errorbar(tau_acf_lf(1:20:end_plot*indc)*1e3,mc(1:20:end_plot*indc)*1e18,Ec(1:20:end_plot*indc)*1e18,'.','MarkerSize',25 ,'LineWidth', 1.5,'Color',colbar, 'DisplayName', 'Experimental ACF');
e.Color = col3;
hold on
plot(tau_acf_lf(1:20:end_plot*indc)*1e3,c0_exp_lf*exp(-tau_acf_lf(1:20:end_plot*indc)/tau0_exp_lf)*1e18, 'LineWidth',3,'Color','k', 'DisplayName',  'Linear fitting')

plot(tau_nl(1:20:end_plot*indc)*1e3,c0_exp_nl*exp(-tau_nl(1:20:end_plot*indc)/tau0_exp_nl)*1e18, '--', 'LineWidth',3,'Color','r', 'DisplayName',  'Non -linear fitting')


ntaus=2.5;
plot([tau0_exp_lf*ntaus,tau0_exp_lf*ntaus]*1e3,[-0.5,400],'--k', 'HandleVisibility','off')
text(tau0_exp_lf*ntaus*partau0*1e3,0.4e2,[num2str(ntaus),'$\tau_{{\rm ot},x}$'],'Interpreter','latex','FontSize',30)
xlim([-0.2 9]);
ylim([-0.1, 3]*1e2)

if aa==5
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5, 'FontSize',25);
else 
   set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5, 'FontSize',25,'Yticklabel',[]);
    end
 
 
 
 xlabel('$\tau(\rm ms)$','Interpreter','Latex', 'FontSize',30)

 
 %% abc
xwi = 400;    % width of the plot square
bx1 = 120;     % extra space at the left
bx2 = 20;     % extra space at the right

Xpix = 3*xwi+bx1+3*bx2;  % total

ywi = 300;    % length riquadro con funzione
by1 = 110;     % extra space below
by2 = 70;     % extra space up

Ypix = 1*by1+1*ywi+1*by2;  % larghezza figura in pixel


 if aa==5
 ylabel('$\rm ACF(    \rm {nm^2})$','Interpreter','Latex','FontSize',30, 'FontName', 'TimesNewRoman')
%  legend ({'a','b','c'},'Box','off','Position',[2*bx1 Ypix/2+20 0 0])
LL= legend ({'Experimental ACF','Linear fitting','Non-linear fitting'},'Box','off','Position',[0.25 0.6 0.1 0.2])

LL.FontSize = 18
 end
 axes('Position',[(0) 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
hold on

xt = [bx1-100,bx1+xwi+bx2,bx1+2*(xwi+bx2)];
yt = [ by1+ywi+30,by1+ywi+30,by1+ywi+30];
str = {'\bf a','\bf b','\bf c'};
text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])

hold off
% 

end
