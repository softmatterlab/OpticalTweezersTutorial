% Auto correlation analysis for Nexp

% Initialization of the workspace
close all

clear;

%close all;

addpath acf

nfile='figs_acf';

figure(1);

set(gcf,'Position',[150 300 800 600])

axes('OuterPosition',[0 0 1 1])

ax1=gca;

subs=1;

kb=1.38064852e-23;

nps=[2,4,6];

col1=[0,0.6,0.8;0.9,0.6,0.2;0.4,0.8,0.4];


% Load data file
for j=1:length(nps)
    
    p=nps(j);
    
    load(['Fig9_1P' num2str(p) '_S.mat'],'datafile','T');
    
    load([nfile filesep datafile '_fit_pot_subs' num2str(subs)  '.mat']);
    
    tau0=(kb*T/(D_acf))/k_acf;
    
    % plots
    axes(ax1)
        
    errorbar(tau(1:20:6*indc),mc(1:20:6*indc)*1e12,Ec(1:20:6*indc)*1e12,'o','Color',col1(1,:),'MarkerSize',5)
    
    hold on
    
    plot(tau(1:6*indc),kb*T/k_acf*exp(-tau(1:6*indc)/tau0)*1e12,'Color',col1(j,:),'LineWidth',1.5)
    
    xlabel('$\tau(s)$','Interpreter','latex')
    
    ylabel('$C_x(\mu \textrm{m}^2)$','Interpreter','latex')
    
    set(gca,'FontSize',20)
    
    axis([0,0.006,-0.1e-4,3.2e-4])
    %
    disp('...')
    
    disp('Autocorrelation function analysis')
    
    disp(['k_acf: ' num2str(k_acf) '+-' num2str(Ek_acf)])
    
    disp(['D_acf: ' num2str(D_acf) '+-' num2str(ED_acf)])
end

savePDF(gcf,[nfile filesep nfile '.pdf'])

