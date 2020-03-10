clear all
clc
close all
figure(1)
DATA=load('data_feeding_cooling_vacuum.mat');
DATAT=struct2table(DATA);

for j=1:size(DATAT,2)
    datai=DATAT.(j);
    plot(datai.freq, datai.psd, 'DisplayName', datai.name, 'LineWidth', 2);
    
    hold on
end
 %legend
 ylabel('psd (pm^2/Hz)')
xlabel('f (kHz)')
set(gca, 'FontSize', 25)
 hold off
 
 
 %%
 
 clear all

figure(2)
DATA=load('data_heterodyne_spectra_vacuum.mat')
DATAT=struct2table(DATA);

for j=1:size(DATAT,2)
    datai=DATAT.(j);
    plot(datai.freq, datai.psd, 'DisplayName', datai.name, 'LineWidth', 2);
    
    hold on
end
 %legend
 ylabel('psd (pm^2/Hz)')
xlabel('f (kHz)')
set(gca, 'FontSize', 25)
 hold off
 %%
 clear all
 DATA=load('data_temperature.mat')
DATAT=struct2table(DATA);
figure()

for j=1:size(DATAT,2)

    datai=DATAT.(j);
    if j==1
        
    loglog(datai.gamma2pi, datai.T,'s',  'DisplayName', datai.name, 'LineWidth', 2);
    hold on
    
    elseif j==4
        
    loglog(datai.gamma2pi, datai.T,'*',  'DisplayName', datai.name, 'LineWidth', 2);
    else
    loglog(datai.gamma2pi, datai.T,  'DisplayName', datai.name, 'LineWidth', 2);
    end

end
legend
hold off
ylabel('T (mK)')
xlabel('gamma/2pi (kHz)')


%%
clear all
DATA=load('data_pressure.mat')

DATAT=struct2table(DATA);
figure(4)

for j=1:size(DATAT,2)
     datai=DATAT.(j);
     semilogy(datai.Time, datai.Pressure, 'DisplayName', datai.name)
     hold on
end

legend

hold off