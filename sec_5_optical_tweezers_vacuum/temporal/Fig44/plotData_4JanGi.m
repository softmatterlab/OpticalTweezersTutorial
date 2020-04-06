clear all;
close all;

fsample = 625000;
fsize = 15;

PSDall = zeros(32769,3);

fprintf('#################\n###########\n remember to enter correct pressure\n ##########\n###############\n');
% gas pressure during measurement in Pascal, which is mbar x 100
Pgas = 1.9E2;

% get calibration to convert voltage to particle position

% first, get all data from calibration measurements
for j=0:19 
    if j<10
        jstring = ['0',num2str(j)];
    else
        jstring = num2str(j);
    end;
    %cd 'TimeTraces'
    filenm = ['17233_p=1p9mbar_V=10V+0.1V_f=45p3kHz_0',jstring];
    %M = csvread([filenm,num2str(j),'.tt']);
    m = (csvread([filenm,'.tt']))';    % read in data
    %cd ..
    [PSDm, absFFTout, FFTm, fm]=MartinFFTnormalized(m, fsample);
    PSDall = PSDall+PSDm;
    var_m(j+1,:) = var(m);
end;

varMean = mean(var_m)   % this variance can be used to gauge amplitudes of signal^2 over time
PSDmean = PSDall/(j+1);



kb = 1.38E-23;  % Boltzmann constant
T = 293;    % room temperature in K
rhoGlass = 2.2E3; % glass density in kg/m^3
eta = 18.27E-6;
dm = 0.372E-9;


% rpart = 132E-9/2;   % particle radius
% Vpart = 4*pi/3*rpart^3;
% mpart = rhoGlass*Vpart;


% fit Lorentzian to data
% select range of data to fit to
fitRangeHz2 = [40,50]*1E3;   % fit range in Hertz
fitRange = round(fitRangeHz2/(fm(2)-fm(1)));

fitRange_excludeHz = [45200,45400];
fitRange_ex = round(fitRange_excludeHz/(fm(2)-fm(1)));

dataToFit = PSDmean(:,1);
dataToFit(fitRange_ex(1)+1:fitRange_ex(2)) = dataToFit(fitRange_ex(1));
dataToFit = dataToFit(fitRange(1)+1:fitRange(2),1);

startingGuess = [1E14 , 1000, 45000, 0];
%[estimates, sse, model] = fitLorentz(fm', PSDmean(:,2), startingGuess);
[estimates, sse, model] = fitLorentzErik((fitRange(1)+1:fitRange(2))'*(fm(2)-fm(1)), dataToFit, startingGuess);
estimates1 = estimates;
[sse, FittedCurve] = model(estimates);

%calib1 = sqrt(kb*T/(pi*estimates(1)*mpart))
freq1 = estimates(3)

% derive the mass from the particle radius using Erik's instructions
% particle radius
r1 = 0.619*9*pi/sqrt(2)*eta*dm^2*Pgas/(rhoGlass*kb*T*2*pi*abs(estimates(2)))
% get derived particle mass
m1 = 4*pi/3*r1^3*rhoGlass;
calibErik1 = 2*pi*sqrt(abs(estimates(1))*pi*m1/(2*abs(estimates(2))*kb*T))
% get calibration using equipartition
calibEquip1 = sqrt(m1*(2*pi*freq1)^2*varMean(1)/kb/T)


% plot for paper

figure;
fm_short = reshape(fm(1:end-1),4,[]);
fm_short = fm_short(1,:);
PSDmean_short = reshape(PSDmean(1:end-1,1),4,[]);
PSDmean_short = PSDmean_short(1,:);
%p=semilogy(fm,PSDmean(:,1)/calibErik1^2,(fitRange(1)+1:fitRange(2))'*(fm(2)-fm(1)),FittedCurve/calibErik1^2);
p=semilogy(fm,PSDmean(:,1)/calibErik1^2,'o',(fitRange(1)+1:fitRange(2))'*(fm(2)-fm(1)),FittedCurve/calibErik1^2);
%p=semilogy(fm_short,PSDmean_short/calibErik1^2,'o');
xlim([40.5 47.5]*1E3);
ylim([2E-19 3E-15]);
set(p,'LineWidth',3);
set(gca,'FontSize',fsize);
xlabel('frequency (Hz)');
ylabel('PSD (m^2/Hz)');
title(['axis 1, Lorentzian fit, d1=',num2str(2*round(r1*1E9)),'nm']);
set(gcf,'PaperPositionMode','auto');
print(gcf,'-depsc',['PSDsAxis1_fit.eps']);




return;
