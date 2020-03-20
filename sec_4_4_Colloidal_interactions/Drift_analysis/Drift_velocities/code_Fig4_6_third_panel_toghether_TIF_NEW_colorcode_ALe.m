%% INITIALIZATION
clear all; close all; clc;


TempeK=300;

calc_2D = true;
%calc_2D = false;

% % fs = 24;
% % fsa = 10;
% % fst = 16;
% % fsl = 16;

bw = 0.02; % um

ta=1.5;
dirtif = ['.' filesep 'tif_xi_list_mod' filesep];
dirtif = ['.' filesep];

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
xi_list = [1 2 4  5  8  9 ];


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
        colo{abs(const_ind_col-3)}=[1    0.7  0   ];
        colo{abs(const_ind_col-4)}=[0.1  0.8    0.05];
        colo{abs(const_ind_col-5)}=[0    0.75 0.9];
        colo{abs(const_ind_col-6)}=[0    0    1   ];
    end
    
    for i = 1:numel(colo)
        colo_light{i} = 0.4*colo{i}+0.6;
    end
    
    for i = 1:numel(colo)
        colo_dark{i} = 0.4*colo{i};
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
filepath=['..' filesep 'LR_d' num2str(trapdist*1e+9) filesep];


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

T_exp_sim = T_exp(1,:);

rbins_theo = rbins_tot;



check_for_low_stat = true;
check_for_err_expe = false;
check_for_err_sim = false;

plot_lowstat = true;  string_kind = '_lowstat';
plot_lowstat = false; string_kind = ''; %'_shift';

plot_diffu_lines = true;


% points to represent
x_range_1_125 = [1 1.25]*(2*R_um);

i0 = find(rbins0==Dcoll_um);
i0_tot = find(rbins_tot==Dcoll_um);

i1 = find(rbins0==1.25*Dcoll_um);
i1_tot = find(rbins_tot==1.25*Dcoll_um);


%% produce TIF figures approrpiately (to avoid rastering)


ywi0 = 160;
xwi0 = 220;

ywi_tif = 5*ywi0;
xwi_tif = 5*xwi0;


numrow = numel(k_list);
N_low_stat_thresh = 2000;



%% figures

numrow = numel(k_list);


lfs = 25;
fslab = 25;
stwdth = 0.6;
hdwdth = 7;
hdlgth = 16;
hdnd = 10;

dyl = 4;
dxl = 4;


%

xs = 2*[35   60  60];
ys = 2*[8   20  10];
%
ywi = ywi0; %80;
xwi = xwi0; %120;

totXpix = sum(xs)+2*xwi;

% totXpix=1400;
totYpix = sum(ys([1 3]))+numrow*ywi+(numrow-1)*ys(2);


xwi1 = ywi;
ywi1 = ywi;

totXpix1 = sum(xs)+2*xwi+xwi1+xs(2);
totXpix1=700;
totYpix1 = 820


hf_drift = figure('Position',[10 10 totXpix1 totYpix1]);



%
s0 = double('a');
% for i=1:numrow
%     lettere3{i,1} = char(s0+(i-1)*3);
%     lettere3{i,2} = char(s0+(i-1)*3+1);
%     lettere3{i,3} = char(s0+(i-1)*3+2);
%     lettere2{i,1} = char(s0+(i-1)*2);
%     lettere2{i,2} = char(s0+(i-1)*2+1);
% end
letfs = 18;
letfs = 18;
let_x0 = -50;
let_y0 = ywi-letfs+6;


for i=[1 3 4 6]
    
%     ah_dri_shade{i,1} = axes('Position',[xs(1) 0 xwi 0]/totXpix1 + [0 ys(1)+(numrow-i)*(ys(2)+ywi) 0 ywi]/totYpix1,'ydir','reverse');
    
    
    if i ==1
    ah_dri_shade{i,1} = axes('Position',[xs(1) 0 xwi 0]/totXpix1 + [0 ys(1)+3*(ys(2)+ywi) 0 ywi]/totYpix1,'ydir','reverse');
    elseif i ==3
        
        
    ah_dri_shade{i,1} = axes('Position',[xs(1) 0 xwi 0]/totXpix1 + [0 ys(1)+2*(ys(2)+ywi) 0 ywi]/totYpix1,'ydir','reverse');
    elseif i ==4
        
             
    ah_dri_shade{i,1} = axes('Position',[xs(1) 0 xwi 0]/totXpix1 + [0 ys(1)+(ys(2)+ywi) 0 ywi]/totYpix1,'ydir','reverse');
            
            elseif i ==6
        
             
    ah_dri_shade{i,1} = axes('Position',[xs(1) 0 xwi 0]/totXpix1 + [0 ys(1) 0 ywi]/totYpix1,'ydir','reverse');
            
            end           
        
        
    % % % %     x_range_driftpara = [2.05 2.65];
    % % % %     y_range_driftpara = [-0.8 0.8]*1.0;
    % % % %     xlim(x_range_driftpara);
    % % % %     ylim(y_range_driftpara);
    x_range_driftpara_tiff = [1 xwi_tif];
    y_range_driftpara_tiff = [1 ywi_tif];
    xlim(x_range_driftpara_tiff);
    ylim(y_range_driftpara_tiff);
    ima = imread([...dirtif 
        'dri_para_' num2str(i) '_tif_new.tif']);
    image(ima)
    axis off; box off; grid off;
    %axis on; box on;
    
%     ah_dri_ax{i,1} = axes('Position',[xs(1) 0 xwi 0]/totXpix1 + [0 ys(1)+(numrow-i)*(ys(2)+ywi) 0 ywi]/totYpix1);
    
    
    if i ==1
    ah_dri_ax{i,1} = axes('Position',[xs(1) 0 xwi 0]/totXpix1 + [0 ys(1)+3*(ys(2)+ywi) 0 ywi]/totYpix1);
    elseif i ==3
        
        
    ah_dri_ax{i,1} = axes('Position',[xs(1) 0 xwi 0]/totXpix1 + [0 ys(1)+2*(ys(2)+ywi) 0 ywi]/totYpix1);
    elseif i ==4
        
             
   ah_dri_ax{i,1} = axes('Position',[xs(1) 0 xwi 0]/totXpix1 + [0 ys(1)+(ys(2)+ywi) 0 ywi]/totYpix1);
            
            elseif i ==6
        
             
   ah_dri_ax{i,1} = axes('Position',[xs(1) 0 xwi 0]/totXpix1 + [0 ys(1) 0 ywi]/totYpix1);
            
            end           
    xlim([0 xwi])
    ylim([0 ywi])
    axis off
    
    
%     ah_dri{i,1} = axes('Position',[xs(1) 0 xwi 0]/totXpix1 + [0 ys(1)+(numrow-i)*(ys(2)+ywi) 0 ywi]/totYpix1);
    
    
    
    
    if i ==1
    ah_dri{i,1} = axes('Position',[xs(1) 0 xwi 0]/totXpix1 + [0 ys(1)+3*(ys(2)+ywi) 0 ywi]/totYpix1);
    elseif i ==3
        
        
    ah_dri{i,1} = axes('Position',[xs(1) 0 xwi 0]/totXpix1 + [0 ys(1)+2*(ys(2)+ywi) 0 ywi]/totYpix1);
    elseif i ==4
        
             
    ah_dri{i,1} = axes('Position',[xs(1) 0 xwi 0]/totXpix1 + [0 ys(1)+(ys(2)+ywi) 0 ywi]/totYpix1);
            
            elseif i ==6
        
             
   ah_dri{i,1} = axes('Position',[xs(1) 0 xwi 0]/totXpix1 + [0 ys(1) 0 ywi]/totYpix1);
            
            end           
    
    
    
    
    x_range_driftpara = [2.05 2.65];
    if i<=4
        y_range_driftpara = [-0.8 0.8]*1.0;
    else
        y_range_driftpara = [-1.4 0.2]*1.0;
    end
    xlim(x_range_driftpara);
    ylim(y_range_driftpara);
    axis off; box off; grid off;
    %axis on; box on;
    
%     ah_dri_shade{i,2} = axes('Position',[xs(1)+xs(2)+xwi 0 xwi 0]/totXpix1 + [0 ys(1)+(numrow-i)*(ys(2)+ywi) 0 ywi]/totYpix1,'ydir','reverse');
    
    
    
    
     if i ==1
     ah_dri_shade{i,2} = axes('Position',[xs(1)+xs(2)+xwi 0 xwi 0]/totXpix1 + [0 ys(1)+3*(ys(2)+ywi) 0 ywi]/totYpix1,'ydir','reverse');
    elseif i ==3
        
        
     ah_dri_shade{i,2} = axes('Position',[xs(1)+xs(2)+xwi 0 xwi 0]/totXpix1 + [0 ys(1)+2*(ys(2)+ywi) 0 ywi]/totYpix1,'ydir','reverse');
    elseif i ==4
        
             
    ah_dri_shade{i,2} = axes('Position',[xs(1)+xs(2)+xwi 0 xwi 0]/totXpix1 + [0 ys(1)+(ys(2)+ywi) 0 ywi]/totYpix1,'ydir','reverse');
            
            elseif i ==6
        
             
    ah_dri_shade{i,2} = axes('Position',[xs(1)+xs(2)+xwi 0 xwi 0]/totXpix1 + [0 ys(1) 0 ywi]/totYpix1,'ydir','reverse');
            
            end           
    
    
    
    
    % % % %     x_range_driftperp = [2.05 2.65];
    % % % %     y_range_driftperp = [-0.8 0.8]*1.0;
    % % % %     xlim(x_range_driftperp);
    % % % %     ylim(y_range_driftperp);
    x_range_driftperp_tiff = [1 xwi_tif];
    y_range_driftperp_tiff = [1 ywi_tif];
    xlim(x_range_driftperp_tiff);
    ylim(y_range_driftperp_tiff);
    ima = imread([dirtif 'dri_perp_' num2str(i) '_tif_new.tif']);
    image(ima)
    axis off; box off; grid off;
    %axis on; box on;
    
%     ah_dri_ax{i,2} = axes('Position',[xs(1)+xs(2)+xwi 0 xwi 0]/totXpix1 + [0 ys(1)+(numrow-i)*(ys(2)+ywi) 0 ywi]/totYpix1);
    
    
    
    
    
     if i ==1
     ah_dri_ax{i,2} = axes('Position',[xs(1)+xs(2)+xwi 0 xwi 0]/totXpix1 + [0 ys(1)+3*(ys(2)+ywi) 0 ywi]/totYpix1);
    elseif i ==3
        
        
    ah_dri_ax{i,2} = axes('Position',[xs(1)+xs(2)+xwi 0 xwi 0]/totXpix1 + [0 ys(1)+2*(ys(2)+ywi) 0 ywi]/totYpix1);
    elseif i ==4
        
             
    ah_dri_ax{i,2} = axes('Position',[xs(1)+xs(2)+xwi 0 xwi 0]/totXpix1 + [0 ys(1)+(ys(2)+ywi) 0 ywi]/totYpix1);
            
            elseif i ==6
        
             
    ah_dri_ax{i,2} = axes('Position',[xs(1)+xs(2)+xwi 0 xwi 0]/totXpix1 + [0 ys(1) 0 ywi]/totYpix1);
            
            end           
    
    
    
    xlim([0 xwi])
    ylim([0 ywi])
    axis off
    
   
    
    
    
    
    
%     ah_dri{i,2} = axes('Position',[xs(1)+xs(2)+xwi 0 xwi 0]/totXpix1 + [0 ys(1)+(numrow-i)*(ys(2)+ywi) 0 ywi]/totYpix1);
    
    
    
    
    
     if i ==1
     ah_dri{i,2} = axes('Position',[xs(1)+xs(2)+xwi 0 xwi 0]/totXpix1 + [0 ys(1)+3*(ys(2)+ywi) 0 ywi]/totYpix1);
    elseif i ==3
        
        
   ah_dri{i,2} = axes('Position',[xs(1)+xs(2)+xwi 0 xwi 0]/totXpix1 + [0 ys(1)+2*(ys(2)+ywi) 0 ywi]/totYpix1);
    elseif i ==4
        
             
    ah_dri{i,2} = axes('Position',[xs(1)+xs(2)+xwi 0 xwi 0]/totXpix1 + [0 ys(1)+(ys(2)+ywi) 0 ywi]/totYpix1);
            
            elseif i ==6
        
             
    ah_dri{i,2} = axes('Position',[xs(1)+xs(2)+xwi 0 xwi 0]/totXpix1 + [0 ys(1) 0 ywi]/totYpix1);
            
            end        
    
    
    x_range_driftperp = [2.05 2.65];
    y_range_driftperp = [-0.8 0.8]*1.0;
    xlim(x_range_driftperp);
    ylim(y_range_driftperp);
    axis off; box off; grid off;
    %axis on; box on;
    
    
    
    
    % for axis arrows
    %
    x0_dri(i) = xwi*(Dcoll_um-x_range_driftpara(1))/(x_range_driftpara(end)-x_range_driftpara(1));
    y0_dri_para(i) = ywi*(0-y_range_driftpara(1))/(y_range_driftpara(end)-y_range_driftpara(1));
    y0_dri_perp(i) = ywi*(0-y_range_driftperp(1))/(y_range_driftperp(end)-y_range_driftperp(1));
    %
    xtival = 1+[0:0.05:0.25];
    xti{i} = xwi*(xtival*Dcoll_um-x_range_driftpara(1))/(x_range_driftpara(end)-x_range_driftpara(1));
    for j=1:numel(xtival)
        if mod(j,2)==0
            xtilab{i,j} = num2str(xtival(j),'%4.2f');
        else
            xtilab{i,j} = num2str(xtival(j),'%3.1f');
        end
    end
    %
    if i<=4
        ytival_dri_para = [0:0.1:0.7];
        ytival_dri_para = [-ytival_dri_para(end:-1:2) ytival_dri_para];
        inde0 = find(ytival_dri_para==0);
        yti_dri_para{i} = ywi*(ytival_dri_para-y_range_driftpara(1))/(y_range_driftpara(end)-y_range_driftpara(1));
        for j=1:numel(ytival_dri_para)
            if mod(j-inde0,5)==0
                ytilab_dri_para{i,j} = num2str(ytival_dri_para(j),'%3.1f');
            else
                ytilab_dri_para{i,j} = ''; %num2str(ytival_dri_para(j),'%3.1f');
            end
        end
        ytilab_dri_para{i,inde0} = '0';
    else
        ytival_dri_para = [-1.4:0.1:0.2];
        %ytival_dri_para = [-ytival_dri_para(end:-1:2) ytival_dri_para];
        inde0 = find(ytival_dri_para==0);
        yti_dri_para{i} = ywi*(ytival_dri_para-y_range_driftpara(1))/(y_range_driftpara(end)-y_range_driftpara(1));
        for j=1:numel(ytival_dri_para)
            if mod(j-inde0,5)==0
                ytilab_dri_para{i,j} = num2str(ytival_dri_para(j),'%3.1f');
            else
                ytilab_dri_para{i,j} = ''; %num2str(ytival_dri_para(j),'%3.1f');
            end
        end
        ytilab_dri_para{i,inde0} = '0';
    end
    %
    ytival_dri_perp = [0:0.1:0.7];
    ytival_dri_perp = [-ytival_dri_perp(end:-1:2) ytival_dri_perp];
    inde0 = find(ytival_dri_perp==0);
    yti_dri_perp{i} = ywi*(ytival_dri_perp-y_range_driftperp(1))/(y_range_driftperp(end)-y_range_driftperp(1));
    for j=1:numel(ytival_dri_perp)
        if mod(j-inde0,5)==0
            ytilab_dri_perp{i,j} = num2str(ytival_dri_perp(j),'%3.1f');
        else
            ytilab_dri_perp{i,j} = '';%num2str(ytival_dri_perp(j),'%3.1f');
        end
    end
    %
    ytilab_dri_perp{i,inde0} = '0';
    
end




%% iterations


N_low_stat_thresh = 2000;

for n = [1  3 4 6] %numel(k_list)
    
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
    
    if check_for_low_stat
        i_highstat_expe = [];
        for i=1:numel(NCc_expe{ik})
            if NCc_expe{ik}(i)>N_low_stat_thresh
                i_highstat_expe = [i_highstat_expe i];
            end
        end
    else
        i_highstat_expe = 1:numel(rbins);
    end
    i_lowstat_expe = setdiff(1:numel(rbins),i_highstat_expe);
    
    index0 = find(Dperp_theo==0);
    Dperp_theo(index0(end)) = 0.401*D1;
    Dperp_theo(index0(1:(end-1))) = NaN;
    Dpar_theo(index0(1:(end-1))) = NaN;
    
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
    
    r_int = x_range_1_125;
    
    igood01 = find(rbins_tot>=r_int(1));
    igood02 = find(rbins_tot<=r_int(2));
    igood_tot0 = intersect(igood01,igood02);  %%%% for _tot quantities
    igood01 = find(rbins>r_int(1));
    igood02 = find(rbins<r_int(2));
    igood0 = intersect(igood01,igood02);      %%%% for normal quantities
    
    
    ifull = intersect(i_highstat_expe,igood0);
    
    igood_tot1 = find(diffupar_err_tot<0.05*D0_spherical_bulk);
    igood_tot = intersect(igood_tot0,igood_tot1);
    
    str_Mac_sim_save{ixi} = num2str(T_exp_sim(ixi)*10);
    
    %%% options
    
    options_diffu_line={'Color',colo{mod(n-1,num_col)+1},'Linestyle','-','Linewidth',2};
    options_drift_line={'Color',colo{mod(n-1,num_col)+1},'Linestyle','-','Linewidth',2};
    %    %
    %     options_diffu_para_expe={0.008,'o','Color',colo{mod(n-1,num_col)+1},'MarkerFaceColor',colo_light{mod(n-1,num_col)+1}};
    %     options_drift_para_expe={0.008,'o','Color',colo{mod(n-1,num_col)+1},'MarkerFaceColor',colo_light{mod(n-1,num_col)+1}};
    %     options_diffu_perp_expe={0.008,'v','Color',colo{mod(n-1,num_col)+1},'MarkerFaceColor',colo_light{mod(n-1,num_col)+1}};
    %     options_drift_perp_expe={0.008,'v','Color',colo{mod(n-1,num_col)+1},'MarkerFaceColor',colo_light{mod(n-1,num_col)+1}};
    %
    options_diffu_para_expe={0.008,'o','Color','k','MarkerFaceColor',colo_light{mod(n-1,num_col)+1},'MarkerSize',3};
    options_drift_para_expe={0.008,'o','Color','k','MarkerFaceColor',colo_light{mod(n-1,num_col)+1},'MarkerSize',3};
    options_diffu_perp_expe={0.008,'o','Color','k','MarkerFaceColor',colo_light{mod(n-1,num_col)+1},'MarkerSize',3};
    options_drift_perp_expe={0.008,'o','Color','k','MarkerFaceColor',colo_light{mod(n-1,num_col)+1},'MarkerSize',3};
    
    arrowproperties = {'StemWidth',stwdth,'HeadWidth',hdwdth,'HeadLength',hdlgth,'HeadNode',hdnd,'color',[1 1 1]*0.5,'EdgeColor','none'};
    
    
    %%% update summary figure %%%%%%%%%%%%%%%%%%%%%%% DRIFT SUMMARY
    figure(hf_drift);
    
    %     set(gcf,'CurrentAxes',ah_dri_shade{n,1}); % parallel
    %     hold on;
    %     % simulation, shade
    %     fill(rbins_tot([igood_tot0 igood_tot0(end:-1:1)]),...
    %         [smooth(dripar_mean_tot(igood_tot0)+dripar_err_tot(igood_tot0)); ...
    %         smooth(dripar_mean_tot(igood_tot0(end:-1:1))-dripar_err_tot(igood_tot0(end:-1:1)))]',...
    %         colo_light{mod(n-1,num_col)+1},'Edgecolor','none');
    %     hold off
    
    set(gcf,'CurrentAxes',ah_dri{n,1}); % parallel
    set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',65,'XMinorTick','off','Tickdir','in','TickLength',[0.02, 0.01]);
    hold on;
    % simulation
    dripar_mean_highlight = 0*rbins_tot + NaN;
    dripar_mean_highlight(igood_tot) = dripar_mean_tot(igood_tot);
    if n<=4
        plot(rbins_tot,dripar_mean_highlight,options_drift_line{:});
    else
        itiny0 = find(rbins_tot>=Dcoll_um);
        itiny1 = find(rbins_tot<=2.21);
        itiny = intersect(itiny0,itiny1);
        plot(rbins_tot(itiny),dripar_mean_highlight(itiny),options_drift_line{:});
    end
    % experiment
    % myerrorbar(B_expe{ik}(ifull)',dripar_av_expe{ik}(ifull)',dripar_std_expe{ik}(ifull)',options_drift_para_expe{:});
    
    errorbar(B_expe{ik}(ifull)',dripar_av_expe{ik}(ifull)',dripar_std_expe{ik}(ifull)', 'o','Color','k','MarkerFaceColor',colo_light{mod(n-1,num_col)+1}, 'MarkerSize',5 ,'LineWidth', 1.5);
   
    %
    hold off;
    drawnow;
    
    % axes
    set(gcf,'CurrentAxes',ah_dri_ax{n,1}); % parallel
    hold on
    % x axis
    arrow2d(x0_dri(n),y0_dri_para(n),xwi,y0_dri_para(n),arrowproperties{:});
    
    for i=2:numel(xti{n})
        % tick
         plot(xti{n}(i) + [0 0],y0_dri_para(n)+[0 -4],'color',[1 1 1]*0.5,'linewidth',ta)
        % tick label
        if (mod(i,2)==1)&&(n>4)
            text(xti{n}(i),y0_dri_para(n)-lfs-6,xtilab{n,i},'fontsize',lfs,'HorizontalAlignment','center','VerticalAlignment','baseline','Interpreter','latex')
        end
        % tick label
        if (mod(i,2)==1)&&(n<=4)&& (n>1)
            text(xti{n}(i),y0_dri_para(n)-lfs-6-30,xtilab{n,i},'fontsize',lfs,'HorizontalAlignment','center','VerticalAlignment','baseline','Interpreter','latex')
        end
        
     if (mod(i,2)==1)&&(n==1)
            text(xti{n}(i),y0_dri_para(n)-lfs-6-20,xtilab{n,i},'fontsize',lfs,'HorizontalAlignment','center','VerticalAlignment','baseline','Interpreter','latex')
        end    
        
    end
    text(xwi,y0_dri_para(n)+dyl,'$r/d$',...
        'Interpreter','latex','fontsize',fslab,...
        'HorizontalAlignment','left','VerticalAlignment','baseline')
    % y axis
    arrow2d(x0_dri(n),0,x0_dri(n),ywi,arrowproperties{:});
    for i=1:numel(yti_dri_para{n})
        if ~strcmpi(ytilab_dri_para{n,i},'')
            % tick
            plot(x0_dri(n)+[0 -4],yti_dri_para{n}(i) + [0 0],'color',[1 1 1]*0.5,'linewidth',ta)
            
            % tick label
            text(x0_dri(n)-3,yti_dri_para{n}(i),ytilab_dri_para{n,i},'fontsize',lfs,'HorizontalAlignment','right','VerticalAlignment','middle','Interpreter','latex')
        else
            plot(x0_dri(n)+[0 -3],yti_dri_para{n}(i) + [0 0],'color',[1 1 1]*0.5,'linewidth',ta)
        end
    end
    text(x0_dri(n)+dxl,ywi+1,'$v_{\parallel} (\mu{\rm m}/{\rm s})$',...
        'Interpreter','latex','fontsize',fslab,...
        'HorizontalAlignment','left','VerticalAlignment','baseline')
    hold off
%     text(let_x0,let_y0,['{\bf ' lettere{n,1} '}'],...
%         'Interpreter','latex','fontsize',letfs,...
%         'HorizontalAlignment','left','VerticalAlignment','baseline')
    drawnow
    %
    
    %     set(gcf,'CurrentAxes',ah_dri_shade{n,2}); % perpendicular
    %     hold on;
    %     % simulation, shade
    %     fill(rbins_tot([igood_tot0 igood_tot0(end:-1:1)]),...
    %         [smooth(driperp_mean_tot(igood_tot0)+driperp_err_tot(igood_tot0)); ...
    %         smooth(driperp_mean_tot(igood_tot0(end:-1:1))-driperp_err_tot(igood_tot0(end:-1:1)))]',...
    %         colo_light{mod(n-1,num_col)+1},'Edgecolor','none');
    %     hold off
    
    set(gcf,'CurrentAxes',ah_dri{n,2}); % perpendicular
    hold on;
    %
    driperp_mean_highlight = 0*rbins_tot + NaN;
    driperp_mean_highlight(igood_tot) = driperp_mean_tot(igood_tot);
    if n<=4
        plot(rbins_tot,driperp_mean_highlight,options_drift_line{:});
    else
        itiny0 = find(rbins_tot>=Dcoll_um);
        itiny1 = find(rbins_tot<=2.21);
        itiny = intersect(itiny0,itiny1);
        plot(rbins_tot(itiny),driperp_mean_highlight(itiny),options_drift_line{:});
    end
    %
%     myerrorbar(B_expe{ik}(ifull),driperp_av_expe{ik}(ifull),driperp_std_expe{ik}(ifull),options_drift_perp_expe{:});
    %
     errorbar(B_expe{ik}(ifull)',driperp_av_expe{ik}(ifull)',driperp_std_expe{ik}(ifull)', 'v','Color','k','MarkerFaceColor',colo_light{mod(n-1,num_col)+1}, 'MarkerSize',5 ,'LineWidth', 1.5);
    hold off;
    
    drawnow;
    
    
    
    
    
    
    
    
    % axes
    set(gcf,'CurrentAxes',ah_dri_ax{n,2}); % parallel
    hold on
    % x axis
    arrow2d(x0_dri(n),y0_dri_perp(n),xwi,y0_dri_perp(n),arrowproperties{:});
    for i=2:numel(xti{n})
        % tick
        
         plot(xti{n}(i) + [0 0],y0_dri_perp(n)+[0 -4],'color',[1 1 1]*0.5,'linewidth',ta)
        
        % tick label
        
               
        if (mod(i,2)==1)&&(n>4)
            text(xti{n}(i),y0_dri_perp(n)-lfs-6,xtilab{n,i},'fontsize',lfs,'HorizontalAlignment','center','VerticalAlignment','baseline','Interpreter','latex')
        end
        
        if (mod(i,2)==1)&&(n<=4)
            text(xti{n}(i),y0_dri_perp(n)-lfs-6-40,xtilab{n,i},'fontsize',lfs,'HorizontalAlignment','center','VerticalAlignment','baseline','Interpreter','latex')
        end
        
        
       
    end
    text(xwi,y0_dri_perp(n)+dyl,'$r/d$',...
        'Interpreter','latex','fontsize',fslab,...
        'HorizontalAlignment','left','VerticalAlignment','baseline')
    % y axis
    arrow2d(x0_dri(n),0,x0_dri(n),ywi,arrowproperties{:});
    for i=1:numel(yti_dri_perp{n})
        if ~strcmpi(ytilab_dri_perp{n,i},'')
            % tick
            plot(x0_dri(n)+[0 -4],yti_dri_perp{n}(i) + [0 0],'color',[1 1 1]*0.5,'linewidth',ta)
            % tick label
            text(x0_dri(n)-3,yti_dri_perp{n}(i),ytilab_dri_perp{n,i},'fontsize',lfs,'HorizontalAlignment','right','VerticalAlignment','middle','Interpreter','latex')
        else
            plot(x0_dri(n)+[0 -3],yti_dri_perp{n}(i) + [0 0],'color',[1 1 1]*0.5,'linewidth',ta)
        end
    end
    text(x0_dri(n)+dxl,ywi+1,'$v_{\perp} (\mu{\rm m}/{\rm s})$',...
        'Interpreter','latex','fontsize',fslab,...
        'HorizontalAlignment','left','VerticalAlignment','baseline')
    hold off
%     text(let_x0,let_y0,['{\bf ' lettere31{n,2} '}'],...
%         'Interpreter','latex','fontsize',letfs,...
%         'HorizontalAlignment','left','VerticalAlignment','baseline')
    
   
    drawnow
    
    
    
end




axes('Position',[(0) 0 totXpix1 0]/totXpix1 + [0 0 0 totYpix1]/totYpix1);  % fa in modo di centrare il riquadro degli assi nella posizione voluta
hold on

xt = [40,160+xwi,40,160+xwi,40,160+xwi,40,160+xwi];
yt = [ totYpix1-20,totYpix1-20,totYpix1-ywi0-20-40,totYpix1-ywi0-20-40,totYpix1-2*(ywi0+40)-20,totYpix1-2*(ywi0+40)-20,totYpix1-3*(ywi0+40)-20,totYpix1-3*(ywi0+40)-20];

 str = {'\bf a','\bf b','\bf c','\bf d','\bf e','\bf f','\bf g','\bf h'};
% str = {'\bf a','\bf b','\bf c'};
text(xt,yt,str,'Interpreter','Latex','FontSize',34)

hold off


axis off

xlim([0 totXpix1])
ylim([0 totYpix1])

% set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'XMinorTick','off','Tickdir','in','TickLength',[0.02, 0.01]);
