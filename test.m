%parpool(4);
%Starting parallel pool (parpool) using the 'local' profile ...
%Connected to the parallel pool (number of workers: 4).
%clear all; close all; tic;

Fr=[1000];
Vis=[0.001];
HydroSize=[60e-9];
H=[10];
Temp=[280,290,300,310,320,330,340];
%Temp=[280,300,320,340];
%Temp=[280];
Core=[15e-9];
Grad=0;  %[-11:0.1:11];

Nnp=10^3;
Ntp=10^5;
Ncyc=10;
Nt=Ncyc*Ntp+1;
Mag=zeros(Nt,3,length(Temp));
%Mag=zeros(Nt-1,3,length(Temp));


Har3=zeros(length(Temp),1); Har5=Har3; Har1=Har3; HarAF1=Har3;
for kT=1:length(Temp);
    tic; 
    kT
    % kFr=1; kV=1; kHyd=1; kH=1; kT=1; kCS=1; kGR=1;
    % calculate the induced magnetization
    %[M,t,dMdt,tt,AField]=BrownSRK4(H,[0,1,0],Fr,Temp(kT),Vis,Nnp,Ncyc,Ntp,HydroSize,Core);
    %[M,t,dMdt,tt,AField]=BrownV2v2old(H,[0,1,0],Fr,Temp(kT),Vis,Nnp,Ncyc,Ntp,HydroSize,Core);
    %[M,t,AField]=Combined_Rotation();
    [M,t,AField]=BrownSRK4(H,[0,1,0],Fr,Temp(kT),Vis,Nnp,Ncyc,Ntp,HydroSize,Core);
    Mag(:,:,kT)=M;
    % find the harmonics
    in=find(AField==max(AField)); inHar=[1,3,5]+0; LC=M(in(end-1)+1:in(end),3); bf=fft(LC);
    Har1(kT)=bf(inHar(1)); Har3(kT)=bf(inHar(2)); Har5(kT)=bf(inHar(3));
    LC=AField(in(end-1)+1:in(end)); bf=fft(LC); HarAF1(kT)=bf(inHar(1));
    toc;
end;

%delete(gcp('nocreate'));

H1h=squeeze(Har1); H3h=squeeze(Har3); H5h=squeeze(Har5); Hafh=squeeze(HarAF1);
Rh=abs(H5h./H3h); Ah=angle(H3h)-angle(H1h);
plot(Temp,Rh,'o:k');