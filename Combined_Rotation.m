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

function [t,M,N,AField]=Combined_Rotation(ft0,ftB,xi0,sig,num,tPts,cycs)

%% DJ
Bv=0.01; %Tesla
f=1000;
[t0,tN,tB,xi0,sig,A,B1,B2]=Combined_parameters(Bv,f)
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

    M(i,:)=mean(m); N(i,:)=mean(n);

    a  = repmat(dot(m,n,2),1,3); %mag easy axis dot prod

    dn = (sig*a.*(m-a.*n)-n/2).*ut + cross(randn(num,3),n).*sqrt(ut); %fixed the sign
    %the first term in previous line is equal to -a*cross(cross(m,n),n), which
    %makes m and n try to align
    n  = n+dn;
    n  = n./repmat(sqrt(n(:,1).^2+n(:,2).^2+n(:,3).^2),1,3);

    xi = XI.*cos(2*pi*i*dt)+2*sig.*a.*n; %total field over time
    %xi = XI+2*sig.*a.*n; %total field over time (static field)

    h  = randn(num,3);
    f1 = cross(xi/al+cross(m,xi),m);
    g1 = cross(h/al+cross(m,h),m);
    mb = m+f1*vt/2+g1*sqrt(vt);

    a2 = repmat(dot(mb,n,2),1,3);
    xb = XI.*cos(2*pi*(i+1)*dt)+2*sig.*a2.*n; %drive field at t+dt
    %xb = XI+2*sig.*a2.*n; %drive field at t+dt

    h2 = randn(num,3);
    f2 = cross(xb/al+cross(mb,xb),mb);
    g2 = cross(h2/al+cross(mb,h2),mb);

    m  = m+(f1+f2)*vt/4+(g1+g2)*sqrt(vt)/2; %average and extra factor of 2 from sde
    %m  = m+(f1)*vt/2+(g1)*sqrt(vt); %average and extra factor of 2 from sde

    m  = m./repmat(sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2),1,3);

    %dm = cross(m,xi/al+cross(m,xi))*vt/2+cross(m,h/al+cross(m,h))*sqrt(vt);
    %m  = m-dm;
    %m  = m./repmat(sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2),1,3);

end


