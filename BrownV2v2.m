%%% DBR 2013 -- Torque balance simulations %%%
%%% DJ edits May 2023

function [M,t,AField,dMdt,tt]=BrownV2v2(Bv,Bs,f,T,visc,N,cycs,tPts,rhy,rco)
%  function [M,t,AField]=BrownV2v2(Bv,Bs,f,T,visc,N,cycs,tPts,rhy,rco)

if nargin<1; Bv = 10; end;        %alternating field in 3rd direction [mT]
if nargin<2; Bs = [0,1,0]; end;   %DJ if nargin<2; Bs = [0,1,0]; end;   %static field in all three directions [mT]
if nargin<3; f = 300; end;       %frequency [Hz]
if nargin<4; T = 340; end;        %temperature [degrees K]
if nargin<5; visc = .001; end;    %viscosity [Pas]
if nargin<6; N = 10^4; end;       %number of particles
if nargin<7; cycs = 5; end; %19; end;      %number of cycles
if nargin<8; tPts =10^4; end;    %time points
if nargin<9;  rhy = 60e-9; end % 60e-9; end;   %hydrodynamic radius [m] (defaults give tE=1.04ms)
if nargin<10; rco = 15e-9; end %15e-9; end;   %core radius [m]

rho  = 3200;   %density [kg/m^3] from data sheet 3.2g/ccm
Bv=Bv/1000;    %change units of fields from mT to T
Bs=Bs/1000;
yp='yes';
Ms   = 76; %70        %saturation magnetization
%Bc   = 0;        %circular field component [mT]

% %% NEW
%  rhy  = 56e-9;   %hydrodynamic radius [m]
%  rco  = 12.5e-9;   %core radius [m]
%  rho  = 5600;    %density [kg/m^3] from data sheet 3.2g/ccm
%  Ms   = 70;      %saturation magnetization [emu/g]
%  gam  = 1.3e9;   %gyromagnetic ratio [Hz/T]
%  al   = 1;       %damping parameter
%  K    = 4000;    %anisotropy constant [J/m^3]
% %%

kT  = (1.38e-23)*T;   %Boltzmann energy
Vhy = (4/3)*pi*rhy^3; %NP hydrodynamic volume [m^3]
Vco = (4/3)*pi*rco^3; %NP core volume [m^3]
%mu  = rho*Ms*Vco;   %magnetic moment [J/T]
%mu=Ms*10^(-18); %See uMod spec sheet
mu=3.2070e-18;
gm  = 6*visc*Vhy;   %drag coefficient [Nsm]
tE  = gm/2/kT;      %Einstein relaxation time [s]
YET = tE/sqrt(1+.21*(mu*Bv/kT)^2);


Dv=kT/gm; Du=mu/gm; %variables
%Dv

%m  = initrand(0,1,N); %init mags randomly
%m   = zeros(N,3); m(:,3) = 1; % init mags in x
m = rand(N,3)*2-1; m = m./repmat(sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2),[1,3]); %size(m)
b = zeros(N,3); b(:,3)=Bv;  % make matrix b in z dir.
%calculate timesteps
if f==0; dt=10^-8; tf=5*tE;
else per=1/f; dt=per/tPts; tf=cycs*per;
end
t=0:dt:tf; AField=zeros(1,length(t));
M   = zeros(length(t),3);
%Du
%Dv

%% Stochastic differential equation loop
in=1;
for i=0:dt:tf
    h=randn(N,3); %stochastic term for fluctuations
    B=b*cos(2*pi*f*i); %drive field over time
    B(:,1)=B(:,1)+Bs(1); B(:,2)=B(:,2)+Bs(2); B(:,3)=B(:,3)+Bs(3);
    AField(in)=B(1,3);

    M(in,:)=sum(m);

    TB=Du*cross(cross(m,B),m); %mag torque
    %Ts=sqrt(2*Dv)*cross(h,m); %stochastic torque
    Ts=sqrt(2*Dv)*[h(:,2).*m(:,3)-h(:,3).*m(:,2),h(:,3).*m(:,1)-h(:,1).*m(:,3),h(:,1).*m(:,2)-h(:,2).*m(:,1)];

    %stochastic diffeq Stratanovich with drift term
    dm = (TB-2*m*Dv)*dt + Ts*sqrt(dt);
    %['Dv', num2str(mean(abs(Ts*sqrt(dt))))]
    % mean(abs((TB-2*m*Dv)*dt))

    %dm = (TB-2*m*Dv)*dt - Ts*sqrt(dt);
    m = m + dm;

    %normalize the magnitude
    nm = sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2);
    m(:,1)=m(:,1)./nm; m(:,2)=m(:,2)./nm; m(:,3)=m(:,3)./nm;

    in=in+1;
end


dt=t(2)-t(1);
tt=t;
%tt=(t(1:end-1)+t(2:end))/2;
dMxdt=(M(2:end,1)-M(1:end-1,1))/dt;
dMydt=(M(2:end,2)-M(1:end-1,2))/dt;
dMzdt=(M(2:end,3)-M(1:end-1,3))/dt;
dMdt(:,1)=dMxdt;
dMdt(:,2)=dMydt;
dMdt(:,3)=dMzdt;

end
