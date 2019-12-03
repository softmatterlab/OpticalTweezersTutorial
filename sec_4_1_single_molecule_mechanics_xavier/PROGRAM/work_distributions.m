%
% Program that computes the free energy of formation of a DNA hairpin.
% 
clear all; close all
%clc

%Adding path to functions
addpath('Functions')

%path to data
dirpath='../RAW_DATA/';
%output path
outpath='../RESULTS/';

%The files contained in dirpath folder are u,f for the unfolding(stretching) and
%folding(releasing) trajectories. They contain 4 columns:
%$1=counter $2=trap distance(nm) $3=time (s) $4=force(pN)
%The files have been sorted in order to have the trap distance order toward
%increasing values.

%number of files in dirpath
nfiles=96;

%KBT in pNnm at Room Temperature =25ºC
KBT=4.114;

%limits for integrating fdx. At l0 all trajectories must be in the folded
%state and at lf they have to be at the unfolded state.
l0=-160;
lf=-60;

%we prepare the vectors for the work computed for each trajectory (the
%integral of fdx
Wu=zeros(1,nfiles);
Wf=zeros(1,nfiles);
%also declare fmax and fmin for each trajectory
fmaxu=zeros(1,nfiles);
fminu=zeros(1,nfiles);
fmaxf=zeros(1,nfiles);
fminf=zeros(1,nfiles);
%We also compute the effective stiffness of the unfolding branch 
keff=zeros(1,nfiles);

for i=1:nfiles
    unfoldfile= char(strcat(dirpath,{'u'},num2str(i),{'.txt'}));
    foldfile= char(strcat(dirpath,{'f'},num2str(i),{'.txt'}));
    
    datau=importdata(unfoldfile);
    dataf=importdata(foldfile);
    %importing force and extension

    xu = datau(:,2);
    fu = datau(:,4);
    %finding indices fullfilling the conditions for integrating
    indexu=find(xu>l0 & xu<lf);
    xuint=xu(indexu);
    fuint=fu(indexu);
    fmaxu(i)=fuint(length(indexu));
    fminu(i)=fuint(1);
    
    xf = dataf(:,2);
    ff = dataf(:,4);
    %finding indices fullfilling the conditions for integrating
    indexf=find(xf>l0 & xf<lf);
    xfint=xf(indexf);
    ffint=ff(indexf);
    fmaxf(i)=ffint(length(indexf));
    fminf(i)=ffint(1);
    
    %to compute the stiffnes we do a linear fit of the first nfit=400points
    nfit=400;
    xfit=xfint(1:nfit);
    ffit=ffint(1:nfit);
    P = polyfit(xfit,ffit,1);
    keff(i)=P(1);
    
    
    %we compute the work for each trajectory and add it to the vector of
    %works
    Wu(i)=trapz(xuint,fuint)/KBT;
    Wf(i)=trapz(xfint,ffint)/KBT;
end
    
%%
% We compute the work histogram for the data analyzed in the previous step
%

%we create the vector of each work distribution. The bining is created by
%using the Rice's rule.

[yyu,xu] = hist(Wu, fix(2*length(Wu)^(1/3)) ); 
[yyf,xf] = hist(Wf, fix(2*length(Wf)^(1/3)) ); 
%normalization: probability distribution functions (pdf)
yu = yyu/trapz(xu,yyu);
yf = yyf/trapz(xf,yyf);

%We adjust the bining so both pdf's have the same bining
W = [min([min(xu),min(xf)]) : 0.5 : max([max(xu),max(xf)]) ]';

%interpolation
Pu = interp1(xu,yu,W);
Pf = interp1(xf,yf,W);

%We find the point at which the distributions coincide (units of pNnm)
Weq=W(find(abs(Pu-Pf) == min(abs(Pu-Pf))));

%We write into a file the pdf of the forward and reverse work:
file_nameF = char(strcat(outpath,'WFpdf.txt'));
dlmwrite(file_nameF, [xf',yf'],'delimiter','\t','precision','%5.5f');
file_nameU = char(strcat(outpath,'WUpdf.txt'));
dlmwrite(file_nameU, [xu',yu'],'delimiter','\t','precision','%5.5f');

%We write into a file the pdf of the interpolated area:
file_namePf = char(strcat(outpath,'Interp_WFpdf.txt'));
dlmwrite(file_namePf, [W,Pf],'delimiter','\t','precision','%5.5f');
file_namePu = char(strcat(outpath,'Interp_Wupdf.txt'));
dlmwrite(file_namePu, [W,Pu],'delimiter','\t','precision','%5.5f');


%%
%Finally we subtract the contributions of the handles and bead, the
%released single-stranded DNA and the hairpin orientation in order to get
%the free-energy of formation of the DNA hairpin


%Optical trap contribution + handles (effective stiffness), also in kBT
Weff=(mean(fmaxf)^2-mean(fminf)^2)/(2*mean(keff))/KBT

%ssDNA contribution using inextensible WLC with elastic parameters:
p=1.35;
l=0.59;
n=44;

WssDNA=IntFdx_WLC(mean(fmaxf),p,l*n,KBT)/KBT

%hairpin orientation assuming d=2nm
d=2.0;
Worient=IntFdx_dipole(d, mean(fmaxf),KBT)/KBT

%Free energy of formation (kBT)
%Values for the sequence of the CD4 molecule are: 51.9KBT (Santalucia, 1998)
G0=Weq-Weff-WssDNA+Worient

%%
%Ploting of Figure 1
%Here we plot  the force distance curves of the unfolding and refolding
%trajectories numbered in plottingfiles
figure(1)
plottingfiles=[1 25 50 75 96]
for k=1:length(plottingfiles)
    i=plottingfiles(k);
    unfoldfile= char(strcat(dirpath,{'u'},num2str(i),{'.txt'}));
    foldfile= char(strcat(dirpath,{'f'},num2str(i),{'.txt'}));
    
    datau=importdata(unfoldfile);
    dataf=importdata(foldfile);
    %importing force and extension
    xu = datau(:,2);
    fu = datau(:,4);
    xf = dataf(:,2);
    ff = dataf(:,4);
    if k==1
        plot(xu,fu,  'r')
        hold on
        plot(xf,ff,  'b')
    else
        plot(xu,fu,  'r')
        plot(xf,ff,  'b')
    end
end
xlabel('\lambda(nm)');
ylabel('f(pN)')
xlim([-180,-40])
pointsx=[l0,lf];
pointsy=[mean(fminu),mean(fmaxf)];
plot(pointsx,pointsy,'o', 'MarkerSize', 15, 'MarkerFaceColor', 'k');
txt1 = '\leftarrow (\lambda_0,f_{min})'; 
text(l0+5,mean(fminu),txt1,'FontSize',14)
txt2 = '(\lambda_1,f_{max}) \rightarrow '; 
text(lf-37,mean(fmaxf),txt2,'FontSize',14)
hold off
%%
figure(2)
%Ploting of Figure 2
%Here we plot the histogram and the probability distribution functions in
%order to obtain figure 2.
histogram(Wf,'DisplayStyle', 'bar', 'Binwidth',1.0,'Normalization', 'pdf','LineWidth',2, 'FaceColor', 'w' ,'EdgeColor','b')
hold on
histogram(Wu,'DisplayStyle', 'bar', 'Binwidth',1.0,'Normalization', 'pdf','LineWidth',2, 'FaceColor', 'w' ,'EdgeColor','r')
plot(W,Pf, '-b','LineWidth',5)
plot(W,Pu, '-r','LineWidth',5)
plot([Weq,Weq],[0,0.3],'--k', 'LineWidth',3)
txt = '\DeltaG_{FU} \rightarrow '; 
text(346,0.27,txt,'FontSize',14)
xlabel('W(k_BT)');
xlim([341,358])
ylabel('P_{U} , P_{F} ');
legend('P_{F}', 'P_{U}');
hold off
%% 


bx1 = 100;     % bordo a sinistra
xwi = 560;    % larghezza riquadro con funzione
bx2 = 30;     % bordino a destra

Xpix = 3*bx1+3*xwi+2*bx2;  % larghezza figura in pixel
Xpix =1400;
by1 = 110;     % bordo in basso
ywi = 500;    % altezza riquadro con funzione
by2 = 50;     % bordo in alto

Ypix = by1+ywi+by2;  % larghezza figura in pixel


col3=[0.00,0.45,0.74];
col4=[0.8500, 0.3250, 0.0980];

figure('Position',[10 20 Xpix Ypix]); % crea la figura
axes('Position',[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
plottingfiles=[1 25 50 75 96]
for k=1:length(plottingfiles)
    i=plottingfiles(k);
    unfoldfile= char(strcat(dirpath,{'u'},num2str(i),{'.txt'}));
    foldfile= char(strcat(dirpath,{'f'},num2str(i),{'.txt'}));
    
    datau=importdata(unfoldfile);
    dataf=importdata(foldfile);
    %importing force and extension
    xu = datau(:,2);
    fu = datau(:,4);
    xf = dataf(:,2);
    ff = dataf(:,4);
    if k==1
        plot(xu,fu,  'r')
        hold on
        plot(xf,ff,  'Color',col3)
    else
        plot(xu,fu,  'r')
        plot(xf,ff,  'Color',col3)
    end
end
xlabel('$\lambda \rm(nm)$','Interpreter','Latex', 'FontSize',30);
ylabel('$f\rm(pN)$','Interpreter','Latex', 'FontSize',30)
xlim([-180,-40])
pointsx=[l0,lf];
pointsy=[mean(fminu),mean(fmaxf)];
plot(pointsx,pointsy,'o', 'MarkerSize', 15, 'MarkerFaceColor', 'k');
plot([-160,-160],[8,11.3825],'--k', 'LineWidth',3)
plot([-60,-60],[8,17.1557],'--k', 'LineWidth',3)
plot([-200,-160],[11.3825,11.3825],'--k', 'LineWidth',3)
plot([-200,-60],[17.1557,17.1557],'--k', 'LineWidth',3)
txt1 = '$(\lambda_0,f_{\rm min})$'; 
text(l0+5,mean(fminu)-0.5,txt1,'Interpreter','latex','FontSize',30)
txt2 = '$(\lambda_1,f_{\rm max})$ '; 
text(lf-60+18,mean(fmaxf)+1,txt2,'Interpreter','latex','FontSize',30)
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01],'XMinorTick','on');

axes('Position',[2*bx1+xwi+50 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix); 

histogram(Wf,'DisplayStyle', 'bar', 'Binwidth',1.0,'Normalization', 'pdf','LineWidth',2, 'FaceColor', 'w' ,'EdgeColor',col3)
hold on
histogram(Wu,'DisplayStyle', 'bar', 'Binwidth',1.0,'Normalization', 'pdf','LineWidth',2, 'FaceColor', 'w' ,'EdgeColor','r')
plot(W,Pf, 'Color',col3,'LineWidth',5)
plot(W,Pu, '-r','LineWidth',5)
plot([Weq,Weq],[0,0.3],'--k', 'LineWidth',3)
txt = '$\Delta G_{\rm FU}$'; 
 
text(345,0.27,txt,'Interpreter','latex','FontSize',30)
xlabel('$W(\rm k_BT)$','Interpreter','latex','FontSize',30);
xlim([341,358])
ylabel('$\rm{pdf}(\it W)$','Interpreter','latex','FontSize',30);
% legend('P_{F}', 'P_{U}','Interpreter','latex','FontSize',30);


LL= legend ({'$\it {W}_{\rm F}$','$\it {W}_{\rm U}$'},'Interpreter','latex','Box','off','Position',[0.87 0.69 0.1 0.2])
LL.FontSize = 18



hold off
set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','on');



axes('Position',[(0) 0 Xpix 0]/Xpix + [0 0 0 Ypix]/Ypix);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
hold on

xt = [bx1-85,bx1+xwi+bx2+8];
yt = [ by1+ywi+30,by1+ywi+30];
str = {'\bf a','\bf b'};
text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 Xpix])
ylim([0 Ypix])
