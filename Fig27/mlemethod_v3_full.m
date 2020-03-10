function  [kx,ky,omegax,omegay,Difx,Dify,Nd,Fx,Fy]=mlemethod_v3_full(Vx,Vy,dt,T,X0,DD,Nlim)
% DD: defines the size of the domain to be computed
Kb=1.38064852e-23;%m^kgs^-2K-1

[lr,N]=size(Vx);
DXo=diff(Vx);
DYo=diff(Vy);

Vx=Vx(1:end-1)-X0(1);
Vy=Vy(1:end-1)-X0(2);

for j=1:N
    Xs=Vx(:,j);
    Ys=Vy(:,j);
    DX=DXo(:,j);
    DY=DYo(:,j);
   % DX=Xs(2:end)-Xs(1:end-1);
   % DY=Ys(2:end)-Ys(1:end-1);
    X=Xs;%(1:end-1);
    Y=Ys;%(1:end-1);
    %filter
    indf=find(sqrt(X.^2+Y.^2)>DD);
    %indf=find((X<-DD/2|X>=DD/2)|(Y<-DD/2|Y>=DD/2));
    %disp(size(indf))
    
    X(indf)=[];
    Y(indf)=[];
    DX(indf)=[];
    DY(indf)=[];
    if length(X)<Nlim
        kx=NaN;
        ky=NaN;
        omegax=NaN;
        omegay=NaN;
        Difx=NaN;
        Dify=NaN;
        Nd=0;
        Fx=NaN;
        Fy=NaN;
    else
    
    Nd=length(X);
    
    ONEM=ones(size(X));
    DXY=1/dt*[DX,DY];
    XY1=[X,Y,ONEM];
    MLE=inv(transpose(XY1)*XY1)*transpose(XY1)*DXY;
    MLE=MLE';
    a=MLE(1,1);
    b=MLE(1,2);
    c=MLE(2,1);
    d=MLE(2,2);
    fx=MLE(1,3);
    fy=MLE(2,3);
    sigmax2=mean((DX/dt-a*X-b*Y-fx).^2);
    sigmay2=mean((DY/dt-c*X-d*Y-fy).^2);
    
    Difx(j)=sigmax2*dt/2;
    Dify(j)=sigmay2*dt/2;
    
    gammax=Kb*T/(sigmax2*dt/2);
    gammay=Kb*T/(sigmay2*dt/2);
    
    
    kx(j)=-gammax*a;
    ky(j)=-gammay*d;
    Fx(j)=gammax*fx;
    Fy(j)=gammay*fy;
    
     J=[MLE(1,1), MLE(1,2); MLE(2,1) , MLE(2,2)];
     MAs=1/2*(J-J');
%      MSim=1/2*(J+J');
%      [V,DD]=eig(MSim);
     omegax(j)=-MAs(1,2);
     omegay(j)=MAs(2,1);
    
     %omegax(j)=b;% revisar los signos
%     omegay(j)=c*gammay;
     end
end

kx=kx';
ky=ky';
Fx=Fx';
Fy=Fy';
omegax=omegax';
omegay=omegay';
Difx=Difx';
Dify=Dify';