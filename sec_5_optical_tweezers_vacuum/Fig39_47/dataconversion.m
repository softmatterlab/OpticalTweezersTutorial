clear all 
close all

load('data.mat') ;
x={xx; yx; zx};
y={xy; yy; zy};
%x conversion
for i=1: 3
  
alp1=-8;
alp2=-2;
xpx=x{i,1};
xpx1=0;
xpx2=260;

alphax=alp1+(xpx-xpx1)/(xpx2-xpx1)*(alp2-alp1);
X{i}=10.^alphax;

%y conversion

alp3=1;
alp4=6;
ypx=y{i,1};
ypx1=0;
ypx2=342;

alphay=alp3+(ypx-ypx1)/(ypx2-ypx1)*(alp4-alp3);
Y{i}=10.^alphay;

figure

scatter (X{i},Y{i})

set(gca,'TickLabelInterpreter','latex',  'yscale', 'log','xscale', 'log')

figure

scatter (xpx,ypx)

end

xx=X{1};
yx=X{2};
zx=X{3};
xy=Y{1};
yy=Y{2};
zy=Y{3};

 save('dataconverted.mat', 'xx', 'xy','yx',  'yy', 'zx', 'zy' )
% set(gca,'TickLabelInterpreter','latex',  'yscale', 'log','xscale', 'log')