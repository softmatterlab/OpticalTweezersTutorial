%% INITIALIZATION
clear all; close all; clc;


set(0, 'defaultFigureRenderer', 'painters')

TempeK=300;

calc_2D = true;
%calc_2D = false;


bw = 0.02; % um


% number of standard deviations for the error
% num_stdev = 2;
num_stdev = 1;


% filepath_expe = ['..' filesep 'Traj16' filesep];
% sottodir_expe='Analysis';
% if ~exist([filepath_expe sottodir_expe filesep],'dir')
%     return
% end


trapdist = 2.38e-6;
eta    = 0.0019;
radius = 1.040e-6;
R_um = radius*1e+6;
Dcoll_um = 2*R_um;
lD  = 12e-9;
lES = 95e-9;


i_e = 5;
i_r = 13;
i_p2m = 9;

%eta=0.002;                %%%%%%%%% CHANGE ACCORDINGLY!!!! cfr line 28
%radius=1.04e-6;           %%%%%%%%% CHANGE ACCORDINGLY!!!! cfr line 29
pix2micron = 0.02800;     %%%%%%%%% CHANGE ACCORDINGLY!!!! cfr line 30

% experimental
load([...filepath_expe sottodir_expe filesep 
    'drift_diffu_dataset' ...
    '_' num2str(i_e) '_' num2str(i_r) '_' num2str(i_p2m) '.mat'],'B',...
    'dri_av','dri_std','diffu_av','diffu_std',...
    'distances','distances_phys','distances_meas','int_TW_list','int_T_list','d_hist',...
    'Cc','Dc','NCc','NDc','Ccpar','Dcpar','Ccperp','Dcperp',...
    'tstep','str_label','T_exp','str_Mac','nt_drift','nt_diffusion','bw','num_files',...
    'D0_spherical_bulk','Xrel','Yrel',...
    'dripar_av','dripar_std','diffupar_av','diffupar_std',...
    'driperp_av','driperp_std','diffuperp_av','diffuperp_std'...
    );

B_expe = B;
dri_av_expe = dri_av;
dri_std_expe = dri_std;
diffu_av_expe = diffu_av;
diffu_std_expe = diffu_std;
distances_expe = distances;
distances_phys_expe = distances_phys;
distances_meas_expe = distances_meas;
int_TW_list_expe = int_TW_list;
int_T_list_expe = int_T_list;
d_hist_expe = d_hist;
Cc_expe = Cc;
Dc_expe = Dc;
NCc_expe = NCc;
NDc_expe = NDc;
Ccpar_expe = Ccpar;
Dcpar_expe = Dcpar;
Ccperp_expe = Ccperp;
Dcperp_expe = Dcperp;
tstep_expe = tstep;
str_label_expe = str_label;
T_exp_expe = T_exp;
str_Mac_expe = str_Mac;
nt_drift_expe = nt_drift;
nt_diffusion_expe = nt_diffusion;
num_files_expe = num_files;
D0_spherical_bulk_expe = D0_spherical_bulk;
dripar_av_expe = dripar_av;
dripar_std_expe = dripar_std;
diffupar_av_expe = diffupar_av;
diffupar_std_expe = diffupar_std;
driperp_av_expe = driperp_av;
driperp_std_expe = driperp_std;
diffuperp_av_expe = diffuperp_av;
diffuperp_std_expe = diffuperp_std;

%% chose T

k_list   = [3 8 10 11 12 13 ];

xi_list = [1 2 4  4  9  9 ];




% asks for index to fit

%% plot options
plot_image = true;

save_images =true;

if plot_image
    col{1}=uint8([1    0    0   ]*255);
    col{2}=uint8([1    0.5  0   ]*255);
    col{3}=uint8([0.8  0.9  0   ]*255);
    col{4}=uint8([0.1  1    0.05]*255);
    col{5}=uint8([0    0.8  0.9 ]*255);
    col{6}=uint8([0    0    1   ]*255);
    col{7}=uint8([0.05 0    0.7 ]*255);
    col{8}=uint8([0.5  0    0.5 ]*255);
    if numel(k_list)>=8
        const_ind_col = 9;
        colo{abs(const_ind_col-1)}=[1    0    0   ];
        colo{abs(const_ind_col-2)}=[1    0.5  0   ];
        colo{abs(const_ind_col-3)}=[1    1    0   ];
        colo{abs(const_ind_col-4)}=[0.1  1    0.05];
        colo{abs(const_ind_col-5)}=[0    0.8  0.9 ];
        colo{abs(const_ind_col-6)}=[0    0    1   ];
        colo{abs(const_ind_col-7)}=[0.05 0    0.7 ];
        colo{abs(const_ind_col-8)}=[0.5  0    0.5 ];
    end
    if numel(k_list)==5
        const_ind_col = 6;
        colo{abs(const_ind_col-1)}=[1    0    0   ];
        colo{abs(const_ind_col-2)}=[1    0.5  0   ];
        colo{abs(const_ind_col-3)}=[1    1    0   ];
        colo{abs(const_ind_col-4)}=[0.1  1    0.05];
        colo{abs(const_ind_col-5)}=[0    0    1   ];
    end
    if numel(k_list)==6
        const_ind_col = 7;
        colo{abs(const_ind_col-1)}=[1    0    0   ];
        colo{abs(const_ind_col-2)}=[1    0.5  0   ];
        colo{abs(const_ind_col-3)}=[1    1    0   ];
        colo{abs(const_ind_col-3)}=[1    1    0   ];
        colo{abs(const_ind_col-4)}=[0.1  1    0.05];
        colo{abs(const_ind_col-5)}=[0    0.75 1];
        colo{abs(const_ind_col-6)}=[0    0    1   ];
    end
    
    for i = 1:numel(colo)
        colo_light{i} = 0.4*colo{i}+0.6;
    end
    
    for i = 1:numel(colo)
        colo_lowstat{i} = 0.4*colo{i}+0.6;
    end
    
    num_col = min([numel(k_list),8]);
end


%pix2micron = 1e+6;    % simulation is in SI units

trapdist = 2.38e-6;
eta = 0.0019; radius = 1.040e-6;
%lD = 13e-9; lES = 95e-9; T_ON = 0.25;
%filepath=['..' filesep 'LR_d' num2str(trapdist*1e+9) filesep];


str_title = 'SIM';

% sottodir='Analysis';
% if ~exist([filepath sottodir filesep],'dir')
%     error('Directory does not exist!!!')
%     return
% end

str_label = {};


% load file
load([...filepath sottodir filesep 
    'histograms.mat'],'d_hist','d_hist0','d_hist1','d_hist_tot',...
    'rbins_hist','rbins0','rbins1','rbins_tot','T_exp');
num_files = size(T_exp,2);

i0 = find(rbins0==Dcoll_um);
i_tot = find(rbins_tot==Dcoll_um);




T_exp_sim = T_exp(1,:);

rbins_theo = rbins_tot;


check_for_low_stat = true;
check_for_err_expe = false;
check_for_err_sim = false;

plot_lowstat = true;  string_kind = '_lowstat';
plot_lowstat = false; string_kind = ''; %'_shift';

plot_diffu_lines = true;



%% figures diffusion

fst = 25;
fsl = 16;
fsa = 10;


we_shift = 0;
vshift_diffu = 0;

totXpix = 430;
totYpix = 300;
% 
% totXpix = 1400;
% totYpix = 330;
%
ys1 = 10;
ys2 = 60;
%
xs1 = 70;
xs2 = 20;
%
ywi = totYpix-ys1-ys2;
xwi = (totXpix-xs1-xs2);

% figure: diffusion parallel
hf_diffu = figure('Position',[10 10 totXpix totYpix]);
%
diffu_ah = axes('Position',[xs1 0 xwi 0]/totXpix + [0 ys2 0 ywi]/totYpix );
x_range_diffu = [1 1.3];
if vshift_diffu==0
    y_range_diffu = [0  1];
    plot_diffu_lines = false;
else
    if vshift_diffu>0
        y_range_diffu = [0 1+numel(k_list)*we_shift*vshift_diffu];
    else
        y_range_diffu = [numel(k_list)*we_shift*vshift_diffu  1];
    end
end
xlim(x_range_diffu);
ylim(y_range_diffu);
axis on; box on; grid off;
xlabel('$r/d$','Interpreter','Latex','fontsize',fsl)
ylabel('$D/D_{\rm 0}$','Interpreter','Latex','fontsize',fsl);
set(diffu_ah,'fontsize',fsa);



%% iterations


N_low_stat_thresh = 2000;

for n = 1:numel(k_list)
    
    ik = k_list(n);
    ixi = xi_list(n);
    
    % experimental parameters chosen
    int_TW_expe = int_TW_list_expe{ik};
    int_T_expe  = int_T_list_expe{ik};
    
    
    % load simulation lines
    load([...filepath sottodir filesep 
        'drift_diffu_mean_err_xi' num2str(T_exp(1,ixi)*10) '.mat'],...
        'rbins','bw','dripar','cdripar','dripar_mean','dripar_err','c_dripar_av',...
        'driperp','cdriperp','driperp_mean','driperp_err','c_driperp_av',...
        'diffupar','cdiffupar','diffupar_mean','diffupar_err','c_diffupar_av',...
        'diffuperp','cdiffuperp','diffuperp_mean','diffuperp_err','c_diffuperp_av',...
        'dripar1','cdripar1','dripar_mean1','dripar_err1','c_dripar_av1',...
        'driperp1','cdriperp1','driperp_mean1','driperp_err1','c_driperp_av1',...
        'diffupar1','cdiffupar1','diffupar_mean1','diffupar_err1','c_diffupar_av1',...
        'diffuperp1','cdiffuperp1','diffuperp_mean1','diffuperp_err1','c_diffuperp_av1',...
        'rbins_theo','Dpar_theo','Dperp_theo','D1','D0_spherical_bulk'...
        );
    
    
    
    Dperp_theo(find(Dperp_theo==0)) = 0.401*D1;
    
%     sottodir_comparison=['Comp_LR_d' num2str(trapdist*1e+6)];
%     if ~exist([filepath_expe sottodir_comparison filesep],'dir')
%         mkdir([filepath_expe sottodir_comparison filesep]);
%     end
    
    
    rbins_tot = reshape([rbins; rbins1],1,2*numel(rbins));
    dripar_mean_tot = reshape([dripar_mean; dripar_mean1],1,2*numel(rbins));
    dripar_err_tot = reshape([dripar_err; dripar_err1],1,2*numel(rbins));
    driperp_mean_tot = reshape([driperp_mean; driperp_mean1],1,2*numel(rbins));
    driperp_err_tot = reshape([driperp_err; driperp_err1],1,2*numel(rbins));
    diffupar_mean_tot = reshape([diffupar_mean; diffupar_mean1],1,2*numel(rbins));
    diffupar_err_tot = reshape([diffupar_err; diffupar_err1],1,2*numel(rbins));
    diffuperp_mean_tot = reshape([diffuperp_mean; diffuperp_mean1],1,2*numel(rbins));
    diffuperp_err_tot = reshape([diffuperp_err; diffuperp_err1],1,2*numel(rbins));
    
    str_legend = num2str(T_exp(1,ixi));
    
    r_int = [2.0 2.8];
    
    igood01 = find(rbins_tot>=r_int(1));
    igood02 = find(rbins_tot<=r_int(2));
    igood_tot0 = intersect(igood01,igood02);  %%%% for _tot quantities
    igood01 = find(rbins>=r_int(1));
    igood02 = find(rbins<=r_int(2));
    igood0 = intersect(igood01,igood02);      %%%% for normal quantities
    
    
    
    
    
    str_Mac_sim_save{ixi} = num2str(T_exp_sim(ixi)*10);
    
    %% visualise histogram
    
    
    
    %%% options
    
    options_diffu_line={'Color',colo{mod(n-1,num_col)+1},'Linestyle','-','Linewidth',2};
    options_drift_line={'Color',colo{mod(n-1,num_col)+1},'Linestyle','-','Linewidth',2};
    options_para_expe={'o','Color','k','MarkerFaceColor',colo{mod(n-1,num_col)+1},'MarkerSize',10,'Linewidth',1};
    options_perp_expe={'v','Color','k','MarkerFaceColor',colo{mod(n-1,num_col)+1},'MarkerSize',10,'Linewidth',1};
%     options_para_expe={'>','Color','k','MarkerFaceColor',colo{mod(n-1,num_col)+1}};
%     options_perp_expe={'d','Color','k','MarkerFaceColor',colo{mod(n-1,num_col)+1}};
    
      
    
    
    %%% update summary figure %%%%%%%%%%%%%%%%%%%%%%% PARALLEL DIFFUSION SUMMARY
    %%% FIGURE
    figure(hf_diffu);
    set(gcf,'CurrentAxes',diffu_ah);
    hold on
    % various cathegorization for EXPE points
    if check_for_low_stat
        i_highstat_expe = [];
        for i=1:numel(NDc_expe{ik})
            if NDc_expe{ik}(i)>N_low_stat_thresh
                i_highstat_expe = [i_highstat_expe i];
            end
        end
    else
        i_highstat_expe = 1:numel(rbins);
    end
    i_lowstat_expe = setdiff(1:numel(rbins),i_highstat_expe);
    if check_for_err_expe
        i_ins_expe = [];
        for i=1:numel(B{ik})
            if (abs(diffupar_av_expe{ik}(i)-diffupar_mean(i))<=2*diffupar_std_expe{ik}(i))
                i_ins_expe = [i_ins_expe i];
            end
        end
    else
        i_ins_expe = 1:numel(rbins);
    end
    i_out_expe = setdiff(1:numel(rbins),i_ins_expe);
    if check_for_err_sim
        i_ins_sim = [];
        for i=1:numel(B{ik})
            if (abs(diffupar_av_expe{ik}(i)-diffupar_mean(i))<=diffupar_err(i))
                i_ins_sim = [i_ins_sim i];
            end
        end
    else
        i_ins_sim = 1:numel(rbins);
    end
    i_out_sim = setdiff(1:numel(rbins),i_ins_sim);
    % various cathegorization for SIM points
    igood_tot1 = find(diffupar_err_tot<0.1*D0_spherical_bulk);
    igood_tot = intersect(igood_tot0,igood_tot1);
    igood1 = find(diffupar_err<0.1*D0_spherical_bulk);
    igood = intersect(igood0,igood1);
    %
    % highlight experiment
    ifull = setdiff(intersect(i_highstat_expe,i_ins_expe),i0);
    ifaded = setdiff(intersect(i_highstat_expe,i_out_expe),i0);
    ifaint = i_lowstat_expe;
    
%     plot(B_expe{ik}(ifull)/Dcoll_um,diffupar_av_expe{ik}(ifull)/2/D1+(n-1)*we_shift*vshift_diffu,options_para_expe{:});
%     plot(B_expe{ik}(ifaded)/Dcoll_um,diffupar_av_expe{ik}(ifaded)/2/D1+(n-1)*we_shift*vshift_diffu,options_para_expe{:});
%     
    
    
    errorbar(B_expe{ik}(ifull)/Dcoll_um,diffupar_av_expe{ik}(ifull)/2/D1+(n-1)*we_shift*vshift_diffu,...
        num_stdev*diffupar_std_expe{ik}(ifull)/2/D1,options_para_expe{:});
    
    errorbar(B_expe{ik}(ifaded)/Dcoll_um,diffupar_av_expe{ik}(ifaded)/2/D1+(n-1)*we_shift*vshift_diffu,...
        num_stdev*diffupar_std_expe{ik}(ifaded)/2/D1,options_para_expe{:});
    
    
    
    if plot_lowstat
        plot(B_expe{ik}(ifaint)/Dcoll_um,diffupar_av_expe{ik}(ifaint)/2/D1+(n-1)*we_shift*vshift_diffu,options_diffu_expe_lowstat{:});
    end
    
    if (n==1) && (we_shift==0)
        plot(rbins_theo/Dcoll_um,2*Dpar_theo/2/D1+(n-1)*we_shift*vshift_diffu,'k','linestyle','-');
    end
    if plot_diffu_lines
        plot(rbins_tot/Dcoll_um,diffupar_mean_tot/2/D1+(n-1)*we_shift*vshift_diffu,'Color',colo{mod(n-1,8)+1});
        diffupar_mean_highlight = 0*rbins_tot + NaN;
        diffupar_mean_highlight(igood_tot) = diffupar_mean_tot(igood_tot);
        plot(rbins_tot/Dcoll_um,diffupar_mean_highlight/2/D1+(n-1)*we_shift*vshift_diffu,options_diffu_line{:});
    end
    hold off;
    drawnow;

    
    figure(hf_diffu);
    set(gcf,'CurrentAxes',diffu_ah);
    hold on;
    % various cathegorization for EXPE points
    if check_for_low_stat
        i_highstat_expe = [];
        for i=1:numel(NDc_expe{ik})
            if NDc_expe{ik}(i)>N_low_stat_thresh
                i_highstat_expe = [i_highstat_expe i];
            end
        end
    else
        i_highstat_expe = 1:numel(rbins);
    end
    i_lowstat_expe = setdiff(1:numel(rbins),i_highstat_expe);
    if check_for_err_expe
        i_ins_expe = [];
        for i=1:numel(B{ik})
            if (abs(diffuperp_av_expe{ik}(i)-diffuperp_mean(i))<=2*diffuperp_std_expe{ik}(i))
                i_ins_expe = [i_ins_expe i];
            end
        end
    else
        i_ins_expe = 1:numel(rbins);
    end
    i_out_expe = setdiff(1:numel(rbins),i_ins_expe);
    if check_for_err_sim
        i_ins_sim = [];
        for i=1:numel(B{ik})
            if (abs(diffuperp_av_expe{ik}(i)-diffuperp_mean(i))<=diffuperp_err(i))
                i_ins_sim = [i_ins_sim i];
            end
        end
    else
        i_ins_sim = 1:numel(rbins);
    end
    i_out_sim = setdiff(1:numel(rbins),i_ins_sim);
    % various cathegorization for SIM points
    igood_tot1 = find(diffupar_err_tot<0.1*D0_spherical_bulk);
    igood_tot = intersect(igood_tot0,igood_tot1);
    igood1 = find(diffuperp_err<0.1*D0_spherical_bulk);
    igood = intersect(igood0,igood1);
    %
    % highlight experiment
    ifull = setdiff(intersect(i_highstat_expe,i_ins_expe),i0);
    ifaded = setdiff(intersect(i_highstat_expe,i_out_expe),i0);
    ifaint = i_lowstat_expe;
    
%     plot(B_expe{ik}(ifull)/Dcoll_um,diffuperp_av_expe{ik}(ifull)/2/D1+(n-1)*we_shift*vshift_diffu,options_perp_expe{:});
%     plot(B_expe{ik}(ifaded)/Dcoll_um,diffuperp_av_expe{ik}(ifaded)/2/D1+(n-1)*we_shift*vshift_diffu,options_perp_expe{:});
    
    errorbar(B_expe{ik}(ifull)/Dcoll_um,diffuperp_av_expe{ik}(ifull)/2/D1+(n-1)*we_shift*vshift_diffu,...
        num_stdev*diffuperp_std_expe{ik}(ifull)/2/D1,options_perp_expe{:});
    
    errorbar(B_expe{ik}(ifaded)/Dcoll_um,diffuperp_av_expe{ik}(ifaded)/2/D1+(n-1)*we_shift*vshift_diffu,...
        num_stdev*diffuperp_std_expe{ik}(ifaded)/2/D1,options_perp_expe{:});
    
    
    
    if plot_lowstat
        plot(B_expe{ik}(ifaint)/Dcoll_um,diffuperp_av_expe{ik}(ifaint)/2/D1+(n-1)*we_shift*vshift_diffu,options_diffu_expe_lowstat{:});
    end
    
    if (n==1) && (we_shift==0)     
        plot(rbins_theo/Dcoll_um,2*Dperp_theo/2/D1+(n-1)*we_shift*vshift_diffu,'k','linestyle','-');
    end
    
    if plot_diffu_lines
        plot(rbins_tot/Dcoll_um,diffuperp_mean_tot/2/D1+(n-1)*we_shift*vshift_diffu,'Color',colo{mod(n-1,8)+1});
        diffuperp_mean_highlight = 0*rbins_tot + NaN;
        diffuperp_mean_highlight(igood_tot) = diffuperp_mean_tot(igood_tot);
        plot(rbins_tot/Dcoll_um,diffuperp_mean_highlight/2/D1+(n-1)*we_shift*vshift_diffu,options_diffu_line{:});
    end
    hold off;
    drawnow;
    
    
    
    
    
    
end

% a few other things

x_shift = -0.2;

figure(hf_diffu);
set(gcf,'CurrentAxes',diffu_ah);
hold on;


text(1.49+x_shift,0.70,'$D_{\perp}$','Interpreter','Latex','fontsize',fst,...
    'HorizontalAlignment','right','VerticalAlignment','bottom');

text(1.49+x_shift,0.30,'$D_{\parallel}$','Interpreter','Latex','fontsize',fst,...
    'HorizontalAlignment','right','VerticalAlignment','top');

hold off;

%set(gca,'TickDir','out');
% set(gca,'fontsize',10);

 set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'XMinorTick','off','Tickdir','in');
drawnow;


figure(hf_diffu);
pix_ah = axes('Position',[xs1 0 xwi 0]/totXpix + [0 ys2 0 ywi]/totYpix );
xlim([0 xwi])
ylim([0 ywi])
hold on;

% two spheres

ra = 12;
thcirc = 0:(pi/20):(2*pi);
xc1U = 3.6*xwi/5;
xc2U = 4.05*xwi/5;
xc1D = 3.6*xwi/5;
xc2D = 4.05*xwi/5;
ycU = 8.7*ywi/10;
ycD = 1.3*ywi/10;

% plot(xc1U+ra*cos(thcirc),ycU+ra*sin(thcirc),'k','LineWidth',2);
% plot(xc2U+ra*cos(thcirc),ycU+ra*sin(thcirc),'k','LineWidth',2);
% plot(xc1D+ra*cos(thcirc),ycD+ra*sin(thcirc),'k','LineWidth',2);
% plot(xc2D+ra*cos(thcirc),ycD+ra*sin(thcirc),'k','LineWidth',2);

fill(xc1U+ra*cos(thcirc),ycU+ra*sin(thcirc),[0 0 0]+0.9,'LineWidth',2);
fill(xc2U+ra*cos(thcirc),ycU+ra*sin(thcirc),[0 0 0]+0.9,'LineWidth',2);
fill(xc1D+ra*cos(thcirc),ycD+ra*sin(thcirc),[0 0 0]+0.9,'LineWidth',2);
fill(xc2D+ra*cos(thcirc),ycD+ra*sin(thcirc),[0 0 0]+0.9,'LineWidth',2);

% plot([xc1U xc2U xc1D xc2D],[ycU ycU ycD ycD],'.k');

% plot([xc1U xc2U],[0 0]+ycU,'k');
% plot([xc1D xc2D],[0 0]+ycD,'k');

sw1 = 0.7;
hw1 = 3;
hl1 = 4;
hn1 = 3;
% la = 7;
% di = 7+2;
la = 8;
di = 4;

arrowproperties2d = {'StemWidth',sw1,'HeadWidth',hw1,'HeadLength',hl1,'HeadNode',hn1};
% arrow2d(xc1U,ycU-ra-di,xc1U,ycU-ra-di+la,'color',[0 0 0]+0.5,arrowproperties2d{:});
% arrow2d(xc2U,ycU+ra+di,xc2U,ycU+ra+di-la,'color',[0 0 0]+0.5,arrowproperties2d{:});
% arrow2d(xc1D-ra-di,ycD,xc1D-ra-di+la,ycD,'color',[0 0 0]+0.5,arrowproperties2d{:});
% arrow2d(xc2D+ra+di,ycD,xc2D+ra+di-la,ycD,'color',[0 0 0]+0.5,arrowproperties2d{:});

arrow2d(xc1U,ycU-ra-di,xc1U,ycU-ra-di-la,'color',[0 0 0]+0.5,arrowproperties2d{:});
arrow2d(xc2U,ycU+ra+di,xc2U,ycU+ra+di+la,'color',[0 0 0]+0.5,arrowproperties2d{:});
arrow2d(xc1D-ra-di,ycD,xc1D-ra-di-la,ycD,'color',[0 0 0]+0.5,arrowproperties2d{:});
arrow2d(xc2D+ra+di,ycD,xc2D+ra+di+la,ycD,'color',[0 0 0]+0.5,arrowproperties2d{:});



hold off;
axis off;


drawnow;



Xpix_new = 700;
Ypix_new = Xpix_new*totYpix/totXpix;

 
 

 set(hf_diffu,'Position',[10 10 Xpix_new Ypix_new]);
  set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'XMinorTick','off','Tickdir','in','TickLength',[0.02, 0.01]);
drawnow


saveas(hf_diffu,['diffu_norm' string_kind '_numstd' num2str(num_stdev)],'fig');




% % figure(hf_diffu);
% % fm = 6;
% % axes('Position',[260 0 fm*16 0]/totXpix+[0 60 0 fm*13]/totYpix);
% % 
% % if exist('inset_640_520.png','file')
% %     im = imread('inset_640_520.png');
% %     image(im)
% % end
% % axis off
% % drawnow

return


