
clear all
close all

load('data.mat');
figure()
scatter (xx,xy);
hold on
scatter (yx,yy);
scatter (zx,zy);



x=xx;
y=xy;
g = fittype('a-b*exp(c*x)');
f0 = fit(x,y,g,'StartPoint',[[1; 50] ; 1]);
xx = linspace(1,200,10);
figure
plot(x,y,'o',xx,f0(xx),'r-');


% x=[1,2,4,6,8]';
% y=[100,140,160,170,175].';
% g = fittype('a-b*exp(-c*x)');
% f0 = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1]);
% xx = linspace(1,8,50);
% figure
% plot(x,y,'o',xx,f0(xx),'r-');