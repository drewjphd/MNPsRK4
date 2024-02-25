%%%  Tests convergence of 1st order and SRK4 magnetization solvers %%%
% DJ Jan 2024

%% Pick an initial time grid size
clear all;
tPts =100;     %time points

%% Set parameter values as desired
Bv = 10;         %alternating field in 3rd direction [mT]
Bs = [0,1,0];    %static field in all three directions [mT]
f = 300;        %frequency [Hz]
T = 300;         %temperature [degrees K]
visc = .001;     %viscosity [Pa-s]
N = 10^4;        %number of particles
cycs = 5;         %19; end;      %number of cycles
rhy = 60e-9;   %hydrodynamic radius [m] (defaults give tE=1.04ms)
rco = 15e-9;  %core radius [m]

%% Test for convergence by iteratively doubling the time grid size and check the "change" 
% iterate until change is less than 0.25%. The change metric is just
% sum of point-wise absolute value Difference divided by sum of pointwise absolute
% value. This should be a robust metric, should not lead to NaNs. Other
% metrics are possible

change=1; %initialize
%initialize
[M,t]=BrownV2v2(Bv,Bs,f,T,visc,N,cycs,tPts,rhy,rco);
figure; plot(t,M); title('1st order');

while change>0.0025;
    tPts=2*tPts
    [Mnext,tnext]=BrownV2v2(Bv,Bs,f,T,visc,N,cycs,tPts,rhy,rco);
    temp=abs(Mnext(1:2:end,:));
    change=mean(sum(abs(Mnext(1:2:end,:)-M))/sum(temp(:)))
    M=Mnext;
    t=tnext;
end
  
figure; plot(t,M); title('1st order');


%Now repeat the experiment for SRK4

change=1; %initialize
tPts=100; %initialize
[M,t]=BrownSRK4(Bv,Bs,f,T,visc,N,cycs,tPts,rhy,rco);
figure; plot(t,M);  title('SRK4');
while change>0.0025;
    tPts=2*tPts
    [Mnext,tnext]=BrownSRK4(Bv,Bs,f,T,visc,N,cycs,tPts,rhy,rco);
    temp=abs(Mnext(1:2:end,:));
    change=mean(sum(abs(Mnext(1:2:end,:)-M))/sum(temp(:)))
    M=Mnext;
    t=tnext;
end
  
figure; plot(t,M); title('SRK4');



