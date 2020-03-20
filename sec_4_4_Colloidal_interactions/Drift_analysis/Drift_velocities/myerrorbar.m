function [ ha ] = myerrorbar( X,Y,E, nleb, varargin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here




ha = errorbar(X,Y,E, varargin{:});
hb = get(ha,'Children'); 

Xdata = get(hb(2),'Xdata');

icentral = 1:9:length(Xdata);
ileft1 = 4:9:length(Xdata);
ileft2 = 7:9:length(Xdata);
iright1 = 5:9:length(Xdata);
iright2 = 8:9:length(Xdata);

% xleft and xright contain the indices of the left and right
%  endpoints of the horizontal lines
xleft = ileft1; xright = iright1;
% modify line length 
Xdata(xleft) = Xdata(icentral) - 0.5*nleb;
Xdata(xright) = Xdata(icentral) + 0.5*nleb;
% xleft and xright contain the indices of the left and right
%  endpoints of the horizontal lines
xleft = ileft2; xright = iright2;
% modify line length 
Xdata(xleft) = Xdata(icentral) - 0.5*nleb;
Xdata(xright) = Xdata(icentral) + 0.5*nleb;


set(hb(2),'Xdata',Xdata);

end

