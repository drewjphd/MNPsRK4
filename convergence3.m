%%%  Tests convergence of 1st order and SRK4 magnetization solvers %%%
% DJ Jan 2024

%% Pick an initial time grid size
clear all;
tPts =10^4;     %time points

%% Set parameter values as desired
Bv = 10;         %alternating field in 3rd direction [mT]
Bs = [0,1,0];    %static field in all three directions [mT]
f = 1000;        %frequency [Hz]
T = 300;         %temperature [degrees K]
visc = .001;     %viscosity [Pa-s]
N = 10^3;        %number of particles
cycs = 5;         %19; end;      %number of cycles
rhy = 60e-9;   %hydrodynamic radius [m] (defaults give tE=1.04ms)
rco = 15e-9;  %core radius [m]

%% Test for convergence by iteratively doubling the time grid size and check the "change" 
% iterate until change is less than 0.25%. The change metric is just
% sum of point-wise absolute value Difference divided by sum of pointwise absolute
% value. This should be a robust metric, should not lead to NaNs. Other
% metrics are possible
tic;
change=1; %initialize
%initialize
%[M,t]=BrownV2v2(Bv,Bs,f,T,visc,N,cycs,tPts,rhy,rco);
%[M,t]=BrownSRK4(Bv,Bs,f,T,visc,N,cycs,tPts,rhy,rco);
[t,M]=Combined_Rotation_SRK4(tPts,N);

figure; plot(t,M); title('1st order');

tPts=linspace(100,10^4,21);
change=zeros(1,21);

for j=1:21
    %tPts=2*tPts
    tic; [tnext,Mnext]=Combined_Rotation(tPts(j),N); toc;
    j

    %temp=abs(Mnext(1:2:end,:));
    %change(j)=mean(sum(abs(Mnext(1:2:end,:)-M))/sum(temp(:)))
    %M=Mnext;
    %t=tnext;
    %M1=spline(t,M(:,1),tnext);
    M2=spline(t,M(:,2),tnext);
    M3=spline(t,M(:,3),tnext);
    count=0;
    for k=1:length(tnext)
        temp=(abs((M3(k)-Mnext(k,3)')/M3(k)));
        if temp<0.5
            change(j)=change(j)+temp; %Skip anomalously large values, cause of division by zero
        else 
            count=count+1;
        end
    end
        change(j)=change(j)/(length(tnext)-count); %Average percentage change per time point
        print('done')
end
  
figure; plot(linspace(100,10^4,21),change,'o--');
toc; 
tic;
%Now repeat the experiment for SRK4

change=zeros(1,21); %initialize
tPts=linspace(100,10^4,21); %initialize
%[M,t]=BrownSRK4(Bv,Bs,f,T,visc,N,cycs,tPts(21),rhy,rco);
%[t,M]=Combined_Rotation_SRK4(tPts,N)

figure; plot(t,M);  title('SRK4');
for j=1:21
    %tPts=2*tPts
    %[Mnext,tnext]=BrownSRK4(Bv,Bs,f,T,visc,N,cycs,tPts(j),rhy,rco);
    tic; [tnext,Mnext]=Combined_Rotation_SRK4(tPts(j),N); toc;
    j

    %M1=spline(t,M(:,1),tnext);
    M2=spline(t,M(:,2),tnext);
    M3=spline(t,M(:,3),tnext);
    count=0;
    for k=1:length(tnext)
        temp=(abs((M3(k)-Mnext(k,3)')/M3(k)));
        if temp<0.5
            change(j)=change(j)+temp; %Skip anomalously large values, cause of division by zero
        else 
            count=count+1;
        end
    end
        change(j)=change(j)/(length(tnext)-count); %Average percentage change per time point
        print('done')
end  

figure; plot(linspace(100,10^4,21),change,'o--');
toc;


