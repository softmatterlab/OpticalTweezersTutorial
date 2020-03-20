function [ fh, ah ] = radial_field_log_quarter_bw_colorcoded( fh, ah, r, Dcoll, forceF, rrange_f, xrange, yrange, rtival, rtilab, stmwth, mf, lfs, F0)
%radial_field Summary of this function goes here
%   Detailed explanation goes here


s2s = r-Dcoll;

%mf = 1e+10;
hdwth = (1/0.4)*stmwth;
hdlgth = (2/0.4)*stmwth;
hdnd = (2/0.4)*stmwth;

figure(fh);
set(fh,'CurrentAxes',ah)
xlim(xrange);
ylim(yrange);
hold on

igood1 = find(s2s<(rrange_f(2)-Dcoll));
igood2 = find(s2s>(rrange_f(1)-Dcoll));
igood = intersect(igood1,igood2);
fmag=forceF(igood);
r2 = s2s(igood);


signF = 0*fmag + 1;
inega = find(fmag<0);
signF(inega) = -1;
% izero = find(forceF==0);
% signF(izero) = 0;


icirc0=find(fmag<0,1,'first');

% % colo_posi=[0.85 0.85 1];
% % colo_nega=[1 0.85 0.85];
colo_posi_arrow=[0 0.5 1];
colo_nega_arrow=[1 0 0];
%colo_axis = [1 1 1]*0.35;

% % if numel(icirc0)
% %     fill(xrange([1 2 2 1 1]),yrange([1 1 2 2 1]),colo_nega,'EdgeColor','none')
% %     fill(r2(icirc0)*cos(0:(pi/50):(2*pi)),r2(icirc0)*sin(0:(pi/50):(2*pi)),...
% %         colo_posi,'EdgeColor','w')
% % else
% %     fill(xrange([1 2 2 1 1]),yrange([1 1 2 2 1]),colo_posi,'EdgeColor','none')
% % end



% now plot the logarithm
fmag=log(abs(forceF(igood)/F0));

inega = find(fmag<0);
fmag(inega) = 0;




theta = (2*pi)*(1:1:24)/24;
phicount=0;
%phi = theta(1:numel((1:10:numel(igood))))/numel((1:10:numel(igood)));
for i=1:20:numel(igood)
    phicount=phicount+1;
    if mod(phicount,2)
        phi(phicount) = 0.5*theta(1); %-0.25*theta(1);
    else
        phi(phicount) = 0.5*theta(1); %-0.25*theta(1);
    end
    for j=1:numel(theta)
    angolo = theta(j)+phi(phicount);
        xin = (r2(i))*cos(angolo);
        yin = (r2(i))*sin(angolo);
        if signF(i)>0
            if mf*fmag(i)>hdlgth
                arrow2d(xin,yin,...
                    xin+mf*fmag(i)*cos(angolo),yin+mf*fmag(i)*sin(angolo),...
                    'stemwidth',stmwth,'headwidth',hdwth,'headlength',hdlgth,'headnode',hdnd,'color',colo_posi_arrow,'EdgeColor','none');
            else
                %plot(xin,yin,'.','color',colo_posi_arrow,'MarkerSize',3);
                arrow2d(xin,yin,...
                    xin+hdlgth*cos(angolo),yin+hdlgth*sin(angolo),...
                    'stemwidth',stmwth,'headwidth',hdwth,'headlength',hdlgth,'headnode',hdnd,'color',colo_posi_arrow,'EdgeColor','none');
            end
        else
            if mf*fmag(i)>hdlgth
                arrow2d(xin,yin,...
                    xin-mf*fmag(i)*cos(angolo),yin-mf*fmag(i)*sin(angolo),...
                    'stemwidth',stmwth,'headwidth',hdwth,'headlength',hdlgth,'headnode',hdnd,'color',colo_nega_arrow,'EdgeColor','none');
            else
                %plot(xin,yin,'.','color',colo_nega_arrow,'MarkerSize',3);
                arrow2d(xin,yin,...
                    xin-hdlgth*cos(angolo),yin-hdlgth*sin(angolo),...
                    'stemwidth',stmwth,'headwidth',hdwth,'headlength',hdlgth,'headnode',hdnd,'color',colo_nega_arrow,'EdgeColor','none');
            end
        end
              
    end
    drawnow
    
    %pause
end

hold off
axis off; box off;

set(gca,'xtick',rtival);
set(gca,'ytick',rtival);
set(gca,'xticklabel',rtilab);
set(gca,'yticklabel',rtilab);
set(gca,'tickdir','out');

set(gca,'fontsize',lfs);

% xlabel('$r-d$ (nm)','Interpreter','Latex','fontsize',lfs);
% ylabel('$r-d$ (nm)','Interpreter','Latex','fontsize',lfs);
% 


drawnow


end

