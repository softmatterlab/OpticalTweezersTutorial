clear all
close all
clc

ipath = pwd;

iname2 = 'TimeTrace_P=50mBar_fdr=135kHz_Vdr=30.2V_PSDdrFit.fig';

ifig2 = fullfile(ipath,iname2);

ipath_Noisefloor = ipath;
iname_NoiseFloor = 'data_NoiseFloor.mat';

load(fullfile(ipath_Noisefloor,iname_NoiseFloor));

fig135 = openfig(ifig2);

axesObjs135    = get(fig135, 'Children');  %axes handles
dataObjs135    = get(axesObjs135, 'Children'); %handles to low-level graphics objects in axes
xdata135    = dataObjs135(2).XData;
ydata135    = dataObjs135(2).YData;
xfit135     = dataObjs135(1).XData;
yfit135     = dataObjs135(1).YData;

figNoise = figure(3)


plotNoise = semilogy(data_NoiseFloor(1,:)/1000,data_NoiseFloor(2,:));

set(fig135, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.3, 0.9*0.35 0.9*0.40]);
set(figNoise, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.1, 0.9*0.35 0.9*0.2]);

xlim = [100 150];
ylim125 = [4 4e3];
ylim135 = [4 2e3];
ylimNoise = [3e-3,1e-2];
xlimInset = [134.8 135.2];
ylimInset = [1e1 2e3];

pause
figure(fig135)
set(gca,'XLim',xlim,'YLim',ylim135,'XTick',[100:10:150],'XGrid','on','YGrid','on');
xlabel('Frequency f (kHz)')
ylabel('Power Spectral Density S_v (bit^2/Hz)')
set(dataObjs135(2),'LineStyle','-','LineWidth',1) %
set(dataObjs135(1),'LineStyle','-','LineWidth',1)
leg2 = legend({'Data','Fit'})
pause
figure(figNoise)
set(gca,'XLim',xlim,'YLim',ylimNoise,'XTick',[100:10:150])
set(gca,'XLim',xlim,'YLim',ylimNoise,'XTick',[100:10:150],'YTick',[3e-3,1e-2],'XGrid','on','YGrid','on')
set(plotNoise,'Color',[000 117 206]/255)
pause
figInset = openfig(ifig2);
axesObjsInset    = get(figInset, 'Children');  %axes handles
dataObjsInset    = get(axesObjsInset, 'Children'); %handles to low-level graphics objects in axes
set(figInset , 'Units', 'Normalized', 'OuterPosition', [0.2, 0.1, 0.15 0.25]);
axis square
set(gca,'XLim',xlimInset,'YLim',ylimInset,'XTick',xlimInset,'YTick',ylimInset)
set(dataObjsInset(2),'Marker','.')
set(dataObjsInset(1),'LineWidth',2)







