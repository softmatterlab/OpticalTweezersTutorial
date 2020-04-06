clear all
close all
col3=[0.00,0.45,0.74];

xwi = 580;    % width of the plot square
bx1 = 200;     % extra space at the left
bx2 = 100;     % extra space at the right
bxm = 20;     % extra space at the midle


Xpix=800;   % total width 
ywi = 400;    % length frame with function

by1 = 100;     % extra space below
by2 = 20;     % extra space up
bym = 30;      %extra space at the midle
Ypix = by1+ywi+by2;  % width in pixel

load('../data/data2_vacuum_0.mat')

% default arguments

kB = 1.38064852e-23; 
radius=136e-9/2;
density=1850;
viscosity = 18.27e-6; % Pascal seconds
c=3e8;
figure('Position',[10 20 Xpix Ypix]);

axes('Position',[bx1 0 xwi 0]/Xpix + [0 by1 0 ywi]/Ypix);  

%establish fittype
fty=fittype('a/(x^2+b)');

%compute guess fit
r0=136e-9/2/2;
mass=particle_mass(r0, density);
[g0, ~] = damping_gas(pressure*100,temperature, radius, density, viscosity);
g0 =g0/(2*pi)
f0=6e4;
fc = f0.^2./g0;
b0=fc^2


fcut=1/4*1/(2*frequency(1));
indf=find(frequency<=fcut)
fwc=frequency(indf);
Pkc=PSD(indf);
a0=Pkc(end)*fwc(end)^2


f=fit(frequency(:),PSD(:), fty,'StartPoint',  [a0, b0])%, 'StartPoint', p0, 'Normalize', 'on','Algorithm', 'Levenberg-Marquardt', 'TolFun', 1e-10, 'TolX', 1e-10)




g0=f0^2/sqrt(f.b)
p(1)=loglog(frequency/1000, PSD/g0, 'LineWidth', 2.5, 'color',col3)
hold on 
p(4)=loglog(frequency/1000,   (f.a./(frequency.^2+f.b))/g0   , '--', 'LineWidth', 2.5, 'Color', 'black')
 
%%

load('../data/data2_vacuum_1.mat')

ftylor=fittype('a./((b-x.^2).^2 +  c*x.^2)');


%guess fit (manual)
g0=1e5;
c0=g0.^2;
f0=6e4;
b0=f0^2;
a0=4.7e19

f=fit(frequency(:),PSD(:), ftylor,'StartPoint',  [a0, b0, c0])%, 'StartPoint', p0, 'Normalize', 'on','Algorithm', 'Levenberg-Marquardt', 'TolFun', 1e-10, 'TolX', 1e-10)

g0=sqrt(f.c);

p(2)=loglog(frequency/1000, PSD/g0, 'LineWidth', 2.5, 'color', [0.9000 0.400 0.1000])
hold on 
loglog(frequency/1000,  (f.a ./((f.b-frequency.^2).^2 +  f.c*frequency.^2))/g0 , '--', 'LineWidth', 2.5, 'Color', 'black')

%%


load('../data/data2_vacuum_2.mat')


% guess fit (manual)
g0=5421.131131169666;
c0=g0.^2;
a0=2.3802449476372086e+18
f0=59009.43175376406;
b0=f0^2;


f=fit(frequency(:),PSD(:), ftylor,'StartPoint',  [a0, b0, c0])
g0=sqrt(f.c);
p(3)=loglog(frequency/1000, PSD/g0, 'LineWidth', 2.5,  'Color', '[0 0.7 0]')
hold on 
loglog(frequency/1000,   (f.a ./((f.b-frequency.^2).^2 +  f.c*frequency.^2))/g0   , '--', 'LineWidth', 2.5, 'Color', 'black')




% plot([x0 x0],[ymin ymax],'--','LineWidth', 2.5, 'Color', 'r')
xlabel('Frequency $f$ (kHz)', 'Interpreter', 'Latex')
ylabel({'Normalized PSD';'$\hat{S}_{vv}(f)/g$ (bit$^2$ {Hz}$^{-2})$'}, 'Interpreter', 'Latex')
 
 LL= legend ([p(1) p(2) p(3) p(4)],'1000 mbar','60 mbar','2.5 mbar','fit','Interpreter','latex','Box','off','Position',[260 300 1 1])%, 'boxoff')

% LL= legend ([p(1) p(4)],'1000 mbar','fit','Interpreter','latex','Box','off','Position',[260 300 1 1])%, 'boxoff')
LL.FontSize = 18;
LL.Box = 'off';
 hold off
 set(gca,'TickLabelInterpreter','latex', 'linewidth',1.5,'FontSize',25,'TickLength',[0.02, 0.01], 'XMinorTick','on','xlim', [1e0,2e2],'ylim', [5e-8,1e-2],...
  'Xtick',[1e0 1e1 1e2 ],'Ytick',[1e-7 1e-5 1e-3 ]);
legend
% saveas(gcf,'Fig35.eps','epsc')
% saveas(gcf,'Fig35.fig')

