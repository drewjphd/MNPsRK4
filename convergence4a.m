%%%  Tests convergence of 1st order and SRK4 magnetization solvers %%%
% DJ Jan 2024

%% Pick an initial time grid size
clear all;
tPts =400;     %time points

%% Set parameter values as desired
Bv = 10;         %alternating field in 3rd direction [mT]
Bs = [0,1,0];    %static field in all three directions [mT]
f = 300;        %frequency [Hz]
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

%change=1; %initialize
%initialize
%[M,t]=BrownV2v2(Bv,Bs,f,T,visc,N,cycs,tPts,rhy,rco);
[M,t]=BrownSRK4(Bv,Bs,f,T,visc,N,cycs,tPts,rhy,rco);
M=M/N;
figure; plot(t,M); title('1st order');
tPts=150;
%tPts=linspace(100,400,21);
change=zeros(1,21);
tic;
NN=linspace(100,10^3,21);
for p=1:21
    N=round(NN(p));
    %for j=1:21
    %tPts=2*tPts
    [Mnext,tnext]=BrownV2v2(Bv,Bs,f,T,visc,N,cycs,tPts,rhy,rco);
    %temp=abs(Mnext(1:2:end,:));
    %change(j)=mean(sum(abs(Mnext(1:2:end,:)-M))/sum(temp(:)))
    %M=Mnext;
    %t=tnext;
    %M1=spline(t,M(:,1),tnext);
    Mnext=Mnext/N;
    M2=spline(t,M(:,2),tnext);
    M3=spline(t,M(:,3),tnext);
    count=0;
    for k=1:length(tnext)
        temp=(abs((M2(k)-Mnext(k,2)')/M2(k)))...
            +(abs((M3(k)-Mnext(k,3)')/M3(k)));
        if temp<0.5
            change(1,p)=change(1,p)+temp/2; %Skip anomalously large values, cause of division by zero
        else
            count=count+1;
        end
    end
    change(1,p)=change(1,p)/(length(tnext)-count); %Average percentage change per time point
    print('done')
    % end
end

figure; plot(NN,change,'o--');
toc;
tic;
%Now repeat the experiment for SRK4

change2=zeros(1,21); %initialize
%tPts=linspace(100,400,21); %initialize
%[M,t]=BrownSRK4(Bv,Bs,f,T,visc,N,cycs,tPts(21),rhy,rco);
figure; plot(t,M);  title('SRK4');


for p=1:21
    N=round(NN(p));
    %    for j=1:21
    %tPts=2*tPts
    [Mnext,tnext]=BrownSRK4(Bv,Bs,f,T,visc,N,cycs,tPts,rhy,rco);
    %M1=spline(t,M(:,1),tnext);
    Mnext=Mnext/N;
    M2=spline(t,M(:,2),tnext)/max(Mnext(:));
    M3=spline(t,M(:,3),tnext)/max(Mnext(:));
    count=0;
    for k=1:length(tnext)
        temp=(abs((M2(k)-Mnext(k,2)')/M2(k)))...
            +(abs((M3(k)-Mnext(k,3)')/M3(k)));
        if temp<0.5
            change2(1,p)=change2(1,p)+temp/2; %Skip anomalously large values, cause of division by zero
        else
            count=count+1;
        end
    end
    change2(1,p)=change2(1,p)/(length(tnext)-count); %Average percentage change2 per time point
    print('done')
    %   end
end

figure; 
plot(NN,change,'x-r','MarkerSize',8,'LineWidth',1);
hold on; plot(NN,change2,'o-b','LineWidth',1);
set(gca,'FontSize',14);
set(gca,'FontWeight','Bold');
title(['Timepoints = ',num2str(tPts)]);
ylim([0 0.14]);
xlabel('Number of NPs, N');
ylabel('Accuracy (mean percentage error)')
legend('SRK4','1st order');
toc;
