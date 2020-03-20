%% INITIALIZATION
clear all; close all; clc;


TempeK=300;
radius=1.03e-6;%+0.5*0.05e-6;
eta=0.003;

calc_2D = true;
%calc_2D = false;


bw = 0.02; % um

% filepath_expe = ['..' filesep 'Traj16' filesep];
% sottodir_expe='Analysis';
% if ~exist([filepath_expe sottodir_expe filesep],'dir')
%     error('EXPERIMENTAL dir does not exist\n');
%     return
% end


trapdist = 2.38e-6;
% %
% filepath_sim=['..' filesep 'd' num2str(trapdist*1e+9) '_D8' filesep];
% sottodir_sim='Analysis';
% if ~exist([filepath_sim sottodir_sim filesep],'dir')
%     error('SIMULATION dir does not exist\n');
%     return
% end

% experimental
% load([filepath_expe sottodir_expe filesep 'drift_diffu_dataset.mat'],'B','dri_av','dri_std','diffu_av','diffu_std',...
%     'distances','distances_phys','distances_meas','int_TW_list','int_T_list','d_hist',...
%     'Cc','Dc','NCc','NDc','tstep','str_label','T_exp','str_Mac','nt_drift','nt_diffusion','num_files','D0_spherical_bulk')

% set these parameters properly
i_e = 5;    % eta = 0.0019;
i_r = 13;   % radius = 1.040e-6;
i_p2m = 9;  % pix2micron = 0.02800;


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



ik_list = [3 10 11 12 13];
%ik_list = [3 10 11 12 13];
const_ind_col = numel(ik_list)+1;

if numel(ik_list)>=8
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
if numel(ik_list)==5
    colo{abs(const_ind_col-1)}=[1    0    0   ];
    colo{abs(const_ind_col-2)}=[1    0.5  0   ];
    colo{abs(const_ind_col-3)}=[1    0.7    0   ];
    colo{abs(const_ind_col-4)}=[0.1  1    0.05];
    colo{abs(const_ind_col-5)}=[0    0    1   ];
end
if numel(ik_list)==6
    colo{abs(const_ind_col-1)}=[1    0    0   ];
    colo{abs(const_ind_col-2)}=[1    0.5  0   ];
    colo{abs(const_ind_col-3)}=[1    0.7    0   ];
    colo{abs(const_ind_col-4)}=[0.1  1    0.05];
    colo{abs(const_ind_col-5)}=[0    0.75 1];
    colo{abs(const_ind_col-6)}=[0    0    1   ];
end

for i = 1:numel(colo)
    colo_light{i} = 0.4*colo{i}+0.6;
end

for i = 1:numel(colo)
    colo_dark{i} = 1*colo{i};
end





load('propagator_traj_choices.mat','i_peak_peak','i_opt_peak','i_peak_opt','i_opt_opt');

for ik = 1:numel(ik_list)
    k = ik_list(ik);
    fh{k}=open(['propagator_' num2str(k) '.fig']);
    ah{k}=get(fh{k},'Children');
end



%% get the info for a new figure




for ik = 1:(numel(ik_list))
    k = ik_list(ik);
    
    posi = get(fh{k},'Position');
    xwi(ik) = posi(3);
    ywi(ik) = posi(4);
    posi = get(ah{k},'Position');
    xwi_norm(ik) = posi(3);
    ywi_norm(ik) = posi(4);
    x_range(ik,:) = get(ah{k},'XLim');
    y_range(ik,:) = get(ah{k},'YLim');
    x_tick{ik} = get(ah{k},'XTick');
    y_tick{ik} = get(ah{k},'YTick');
    x_tick_lab{ik} = get(ah{k},'XTickLabel');
    y_tick_lab{ik} = get(ah{k},'YTickLabel');
end

%% new figure
xwi_ff = xwi.*xwi_norm;
ywi_ff = ywi.*ywi_norm;


% all in a row

dyt = 10;
dyb = 100;
dxl = 100;
dxm = 10;
dxr = 10;

Ypix = dyt + ywi_ff(1) + dyb;
Xpix = dxl + sum(xwi_ff(1:(end-1))) + (numel(ik_list)-2)*dxm + dxr;

ffh = figure('Position',[10 10 Xpix Ypix]);

fsa = 10;
fsl = 16;

for ik = 1:(numel(ik_list)-1)
    
    if ik<(numel(ik_list)-1)
        
        
        k = ik_list(ik);
        
        fah{ik} = axes('Position',[dxl+sum(xwi_ff(1:(ik-1)))+(ik-1)*dxm 0 xwi_ff(ik) 0]/Xpix+...
            [0 dyb 0 ywi_ff(ik)]/Ypix);
        
        kids = get(ah{k},'Children');
        
        for i = numel(kids):-1:1
            copyobj(kids(i),fah{ik});
        end
        
        set(fah{ik},'XLim',x_range(ik,:));
        set(fah{ik},'YLim',y_range(ik,:));
        set(fah{ik},'fontsize',fsa);
        set(fah{ik},'Tickdir','out');
        set(fah{ik},'XMinorTick','on');
        set(fah{ik},'YMinorTick','on');
        set(fah{ik},'XTick',x_tick{ik});
        set(fah{ik},'XTickLabel',x_tick_lab{ik});
        xlabel('$t$ [ms]','Interpreter','Latex','fontsize',fsl);
        if ik == 1
            set(fah{ik},'YTick',y_tick{ik});
            set(fah{ik},'YTickLabel',y_tick_lab{ik});
            ylabel('$r$ [$\mu$m]','Interpreter','Latex','fontsize',fsl);
        else
            set(fah{ik},'YTick',[]);
        end
        
        
        
    else
        k = ik_list(ik+1);
        
        fah{ik} = axes('Position',[dxl+sum(xwi_ff(1:(ik-1)))+(ik-1)*dxm 0 xwi_ff(ik) 0]/Xpix+...
            [0 dyb 0 ywi_ff(ik)]/Ypix);
        
        kids = get(ah{k},'Children');
        
        for i = numel(kids):-1:1
            copyobj(kids(i),fah{ik});
        end
        
        set(fah{ik},'XLim',x_range(ik,:));
        set(fah{ik},'YLim',y_range(ik,:));
        set(fah{ik},'fontsize',fsa);
        set(fah{ik},'Tickdir','out');
        set(fah{ik},'XMinorTick','on');
        set(fah{ik},'YMinorTick','on');
        set(fah{ik},'XTick',x_tick{ik});
        set(fah{ik},'XTickLabel',x_tick_lab{ik});
        xlabel('$t$ [ms]','Interpreter','Latex','fontsize',fsl);
        if ik == 1
            set(fah{ik},'YTick',y_tick{ik});
            set(fah{ik},'YTickLabel',y_tick_lab{ik});
            ylabel('$r$ [$\mu$m]','Interpreter','Latex','fontsize',fsl);
        else
            set(fah{ik},'YTick',[]);
        end
    end
    drawnow;
    
end



Xpix_new = 600;
Ypix_new = Ypix/Xpix*Xpix_new;

fs = 24;


set(ffh,'Position',[10 10 Xpix_new Ypix_new]);

drawnow;






stringalettera{1} = '{\bf a}';
stringalettera{2} = '{\bf b}';
stringalettera{3} = '{\bf c}';
stringalettera{4} = '{\bf d}';
stringalettera{5} = '{\bf e}';
axes('Position',[0 0 1 1]);
hold on
for ik=1:(numel(ik_list)-1)
    text((dxl+sum(xwi_ff(1:(ik-1)))+(ik-1)*dxm+15)/Xpix,...
        (dyb+ywi_ff(ik)-40)/Ypix,stringalettera{ik},'Interpreter','Latex','fontsize',fs,...
        'HorizontalAlignment','left','VerticalAlignment','baseline')
end
hold off
axis off




