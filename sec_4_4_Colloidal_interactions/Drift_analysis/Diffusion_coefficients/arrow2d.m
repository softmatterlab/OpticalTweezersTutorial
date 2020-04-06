function h = arrow2d(x1,y1,x2,y2,varargin)
% ARROW2D   Draw a line with an arrowhead in 2D
%
% H = ARROW2D(X1,Y1,X2,Y2) draws arrow(s) from (X1,Y1) to (X2,Y2)
%   and returns the graphics handle of the arrow(s).
%   X1,Y1,X2,Y2 should be matrices with the same dimensions.
%
% H = ARROW2D(X1,Y1,X2,Y2,'PropertyName',PropertyValue) sets the property
%   PropertyName to PropertyValue. All standard plot properties
%   can be used and also the ARROW2D properties listed below.
%
% ARROW2D properties:
%       Color       Arrow color [default = 'k']
%       StemWidth   Arrow stem width [default = .1]
%       HeadLength  Arrow head length [default = 1]
%       HeadWidth   Arrow head width [default = 1]
%       HeadNode    Arrow head intersection with the arrow stem [default = .5]
%
% See also ARROW3D, PLOT.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

% Permission is granted to distribute ARROW3D with the toolboxes for the book
% "Optical Tweezers", by P. H. Jones, O. M. Marago & G. Volpe 
% (Cambridge University Press, 2015).

% Arrow color
color = 'k';
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'color')
        color = varargin{n+1};
    end
end

% Arrow stem width
swidth = .1;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'stemwidth')
        swidth = varargin{n+1};
    end
end

% Arrow head length
hlength = 1;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'headlength')
        hlength = varargin{n+1};
    end
end

% Arrow head width
hwidth = 1;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'headwidth')
        hwidth = varargin{n+1};
    end
end

% Arrow head intersection with the arrow stem
hnode = .5;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'headnode')
        hnode = varargin{n+1};
    end
end

% Calculates the coordinates of the arrow in a standard reference frame
x1 = reshape(x1,1,numel(x1));
x2 = reshape(x2,1,numel(x2));
y1 = reshape(y1,1,numel(y1));
y2 = reshape(y2,1,numel(y2));

slength = sqrt((x2-x1).^2+(y2-y1).^2); % length

v = ones(1,numel(x1));
x =  [ ...
    0*v; ...
    0*v; ...
    (slength-hnode).*v; ...
    (slength-hlength).*v; ...
    slength.*v; ...
    (slength-hlength).*v; ...
    (slength-hnode).*v; ...
    0*v; ...    
    ];
y = [ ...
    swidth.*v; ...
    -swidth.*v; ...
    -swidth.*v; ...
    -hwidth.*v; ...
    0*v; ...
    hwidth.*v; ...
    swidth.*v; ...
    swidth.*v ...
    ];

% Rotates and translates the arrow
angle = atan2(y2-y1,x2-x1);
cosa = ones(size(x,1),1)*cos(angle);
sina = ones(size(x,1),1)*sin(angle);
xr = x.*cosa - y.*sina;
yr = x.*sina + y.*cosa;

xrt = xr+ones(size(x,1),1)*x1;
yrt = yr+ones(size(x,1),1)*y1;

% Plots the arrow(s)
h = fill(xrt,yrt,color,'EdgeColor',color);

% Sets other properties
for n = 1:2:length(varargin)
    if ~strcmpi(varargin{n},'color') ...
            & ~strcmpi(varargin{n},'stemwidth') ...
            & ~strcmpi(varargin{n},'headlength') ...
            & ~strcmpi(varargin{n},'headwidth') ...
            & ~strcmpi(varargin{n},'headnode')
        set(h,varargin{n},varargin{n+1});
    end
end
