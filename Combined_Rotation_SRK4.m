%%%% Solves combined rotation stochastic eqs %%%%
%%%% DBR Feb 2015 %%%%
%%% 2018 fixed bug in line 40, which misses a negative sign

%inputs:
% frequency of oscillating field [Hz]
% Neel event time [s]
% Brownian relaxation time [s]
% Unitless applied field magnitude []
% Unitless anisotropy field magnitude []
% Number of particles []
% Time points per cycle []
% Number of cycles []
% (ft0=f*t0, ftB=f*tB, and other parameters can be obtained using
% Combined_parameters.m)

%outputs:
% Time points
% The magnetization direction at every time
% The easy axis direction at every time

function [t,M,N,AField]=Combined_Rotation_SRK4(ft0,ftB,xi0,sig,num,tPts,cycs);

%% DJ
Bv=0.01; %Tesla
f=1000; %Hz
[t0,tN,tB,xi0,sig,A,B1,B2]=Combined_parameters(Bv,f); %Generate parameters
ft0=f*t0;
ftB=f*tB;
num=10^4;
tPts=10^5;
cycs=5;
%%

al = 1;    %LLG alpha
m  = initrand(num,0);   %initial magnetizations
%m = repmat([0,0,1],num,1); %all in z direction
n  = m;               %initial easy axes along magnetization

XI = repmat([0 0 xi0],num,1); %unitless field vector

t=linspace(0,cycs,tPts*cycs); dt=t(2); %time points

ut=dt/ftB; vt=dt/ft0; %unitless timesteps (Yipeng: dt/f is the actual time for each step and has unit [s])

M=zeros(length(t),3); N=M; %init mean mag and easy axes

for i=1:length(t)
    AField(i)=cos(2*pi*i*dt);
    tt=i*dt;
    M(i,:)=mean(m); N(i,:)=mean(n);

    a  = repmat(dot(m,n,2),1,3); %mag easy axis dot prod

    % DJ dn = (sig*a.*(m-a.*n)-n/2).*ut + cross(randn(num,3),n).*sqrt(ut); %fixed the sign

    % Compute stochastic increments
    %dW = sqrt(ut) * randn(N,3);
    h  = randn(num,3);
    rr  = randn(num,3);

    nm=zeros(length(n),6); %DJ: A "super-vector" containing both n and m
    nm(:,1:3)=n;
    nm(:,4:6)=m;
    uvt=zeros(1,6); %DJ: Same idea for time. Another super-vector
    uvt(1:3)=ut;
    uvt(4:6)=vt;
    dW=zeros(num,6);
    dW(:,1:3)=sqrt(ut);
    dW(:,4:6)=sqrt(vt);
    %h  = randn(num,3);
    %h = rr;

    % Compute RK4 increments
    k1 = c(tt, nm, sig, a, XI, al) .* uvt + d(nm, h, al, rr) .* dW;
    k2 = c(tt + 0.5 * dt, nm + 0.5 * k1, sig, a, XI, al) .* uvt + d(nm + 0.5 * k1, h, al,  rr).* dW;
    k3 = c(tt + 0.5 * dt, nm + 0.5 * k2, sig, a, XI, al) .* uvt + d(nm + 0.5 * k2, h, al,  rr) .* dW;
    k4 = c(tt + dt, nm + k3, sig, a, XI, al) .* uvt + d(nm + k3, h, al,  rr).* dW;

    % Update n using RK4 formula
    dnm = (1/6) * (k1 + 2 * k2 + 2 * k3 + k4);
    nm=nm+dnm;

    %the first term in previous line is equal to -a*cross(cross(m,n),n), which
    %makes m and n try to align
    n  = nm(:,1:3);
    m  = nm(:,4:6);
    n  = n./repmat(sqrt(n(:,1).^2+n(:,2).^2+n(:,3).^2),1,3);
    m  = m./repmat(sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2),1,3);
end
end

% "helper" functions
function val=c(t, nm, sig, a, XI, al)
val=zeros(length(nm),6);
n=nm(:,1:3);
m=nm(:,4:6);
val(:,1:3)=sig*a.*(m-a.*n)-n/2;
xi = XI.*cos(2*pi*t)+2*sig.*a.*n; %total field over time
f1 = cross(xi/al+cross(m,xi),m);
val(:,4:6) = m+f1/2;
end

function val=d(nm, h, al, rr)
val=zeros(length(nm),6);
n=nm(:,1:3);
m=nm(:,4:6);
val(:,1:3)=cross(rr,n);
val(:,4:6)=cross(h/al+cross(m,h),m);%.*rrr;
end
