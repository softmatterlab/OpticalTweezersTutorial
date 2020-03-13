%% Program to obtain separated trajectories and stiffness of the cell from pulling protocol
% Given a txt file, the program ask to select the number of columns of the file, the columns of Y force, 
% trap A Y distance, trap B Y distance, value of trap stiffness, maximum force, minimum force and number of 
% stiffness intervals,it returns the separated trajectories in diferent files named
% 'filename_#trajectoryT_U.txt'(ex:RBC1_1T_U.txt') if we are going from low forces to high forces and
% 'filename_#trajectoryT_D.txt'(ex:RBC1_1T_D.txt') if we are going from high forces to low forces.
% In these files the first column will be celullar extension and the second column will be Y force. 
% All the trajeectories will be forced to start to (0,0) in order to correct any drift. 
% The second file named 'filename_Stiffness.txt' (ex: RBC1_Stiffness.txt) where first column will
% correspond to the mean Y force and the second one to the Stiffness value.

clear all;close all;clc;

[file,path] = uigetfile('*.txt');
cd (path);
Ncol=input('Please introduce the number of columns :','s');

NColF='%f\t';
for m=1:str2num(Ncol)-2
    NColF=char(strcat(NColF,'%f\t'));
end

NColF=char(strcat(NColF,' %f\n'));
FC=input('Column of Y force :','s');
ADC=input('Column of trap A Y distance :','s');
BDC=input('Column of trap B Y distance :','s');
fid=fopen(file,'r');
data = textscan(fid,NColF, 'CommentStyle', '#');

F=data{str2num(FC)};
F=-F; %Just in case you have it in negative values as example data
AD=data{str2num(ADC)};
BD=data{str2num(BDC)};
D=(AD+BD)/2; % The distance that the trap is moved is an average of the displacement of both lasers (double trap)
kkt=input('Trap stiffness value in pN/nm :','s');
kt=str2num(kkt); %pN/nm trap stiffness
E=D-F/kt;

SS=0; %Artificial status. We assume that we start at minimum forces and we move to higher forces
Ntraj=0;
MaF=input('Maximum force :','s');
MiF=input('Minimum force :','s'); % You have to reach the minimum force each cycle. Check it!
MaxF=str2num(MaF);
MinF=str2num(MiF);
for i=1:length(F)
    if (SS==0) 
        Ntraj=Ntraj+1;
        %File where we will write Y distance and Y force of  UP trajectory
        fid1=fopen(char(strcat('../results/', file(1:end-4),'_',num2str(Ntraj),'T_U.txt')),'a');
        fprintf(fid1,'%s\t%s\n','#Y dist(nm)', 'Y force(pN)');
        E=E-E(i);% We put D=0nm at minimum force
        SS=1;
        fprintf(fid1,'%f\t%f\n',E(i),F(i));
    elseif (SS==1&&F(i)<MaxF) 
        fprintf(fid1,'%f\t%f\n',E(i),F(i));
    elseif (SS==1&&F(i)>MaxF)
        fclose(fid1);
        %File where we will write Y distance and Y force of DOWN trajectory
        fid2=fopen(char(strcat('../results/', file(1:end-4),'_',num2str(Ntraj),'T_D.txt')),'a');
        fprintf(fid2,'%s\t%s\n','#Y dist(nm)', 'Y force(pN)');
        SS=2;
        fprintf(fid2,'%f\t%f\n',E(i),F(i));
    elseif (SS==2&&F(i)>MinF) 
        fprintf(fid2,'%f\t%f\n',E(i),F(i));
    elseif (SS==2&&F(i)<MinF)
        SS=0;
        fclose(fid2);
    end   
end

%% Second part of the program: read trajectories cuted files and extract the stiffness 

% We will not consider the first neither the last trajectory
F=0;
E=0;
P=input('Number of intervals for computing stiffness :','s');
points=str2num(P);

for i=2:Ntraj-1
    %File where we will write Y distance and Y force of  UP trajectory
    fid4=fopen(char(strcat('../results/', file(1:end-4),'_',num2str(i),'T_U.txt')),'r');
    fid5=fopen(char(strcat('../results/', file(1:end-4),'_',num2str(i),'T_D.txt')),'r');
    data2 = textscan(fid4,'%f\t%f\n', 'CommentStyle', '#');
    data3 = textscan(fid5,'%f\t%f\n', 'CommentStyle', '#');
    fclose(fid4);
    fclose(fid5);
    E=data2{1};
    F=data2{2};
    E3=data3{1};
    F3=data3{2};
    ind=round(linspace(MinF,MaxF,points))
    for m=1:length(ind)-1
        LinF=fittype('A*x+B');
        startPoints=[0.001,1];
        index=find((F>ind(m)).*(F<ind(m+1))==1);
        index3=find((F3>ind(m)).*(F3<ind(m+1))==1);
        fLinF=fit(E(index),F(index),LinF,'Start', startPoints)
        Cof=coeffvalues(fLinF);
        CI = confint(fLinF,.95);
        k(i-1,m)=Cof(1);
        MF(i-1,m)=mean(F(index));
        Cof=0;
        CI=0;
        index=0;
        fLinF3=fit(E3(index3),F3(index3),LinF,'Start', startPoints)
        Cof3=coeffvalues(fLinF3);
        CI3 = confint(fLinF3,.95);
        k3(i-1,m)=Cof3(1);
        MF3(i-1,m)=mean(F3(index3));
        Cof3=0;
        CI3=0;
        index3=0;
    end
    F=0;
    E=0;
    F3=0;
    E3=0;
end

fid6=fopen(char(strcat('../results/', file(1:end-4),'_Stiffness_U.txt')),'a');
fid7=fopen(char(strcat('../results/', file(1:end-4),'_Stiffness_D.txt')),'a');
fprintf(fid6,'%s\t%s\t%s\t%s\n','#Y force(pN)', 'Err force(pN)', 'Stiffness(pN/nm)','Err Stiffness(pN/nm)');
fprintf(fid7,'%s\t%s\t%s\t%s\n','#Y force(pN)', 'Err force(pN)', 'Stiffness(pN/nm)','Err Stiffness(pN/nm)');
for i=1:length(k(1,:))
    fprintf(fid6,'%f\t%f\t%f\t%f\n',mean(MF(:,i)),std(MF(:,i)),mean(k(:,i)),std(k(:,i)));
    fprintf(fid7,'%f\t%f\t%f\t%f\n',mean(MF3(:,i)),std(MF3(:,i)),mean(k3(:,i)),std(k3(:,i)));
end
fclose(fid6);
fclose(fid7);
