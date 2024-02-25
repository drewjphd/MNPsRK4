%%% DBR 2013 -- Torque balance simulations %%%

%%% DJ 2023: 
% Upgraded numerical integration method from simple Euler-Maryuma 
% to 4th Order Stochastic Runge-Kutta to improve solution accuracy, stability and speed. 
% Smaller N and larger dt suffices. Formula from Appendix of following paper:
% Hansen, J. A., and C. Penland, 2006:
% Efficient Approximate Techniques for Integrating Stochastic Differential Equations.
% Mon. Wea. Rev., 134, 3006–3014, https://doi.org/10.1175/MWR3192.1.

function [M,t,AField]=BrownSRK4(Bv,Bs,f,T,visc,N,cycs,tPts,rhy,rco)

if nargin<1; Bv = 10; end;        %alternating field in 3rd direction [mT]
if nargin<2; Bs = [0,1,0]; end;   %static field in all three directions [mT]
if nargin<3; f = 300; end;       %frequency [Hz]
if nargin<4; T = 340; end;        %temperature [degrees K]
if nargin<5; visc = 0.001; end;    %viscosity [Pas]
if nargin<6; N = 10^4; end;       %number of particles
if nargin<7; cycs = 5; end;        %number of cycles
if nargin<8; tPts =10^3; end;    %time points
if nargin<9;  rhy = 60e-9; end % 60e-9; end;   %hydrodynamic radius [m] (defaults give tE=1.04ms)
if nargin<10; rco = 15e-9; end %15e-9; end;   %core radius [m]

rho  = 3200;   %density [kg/m^3] from data sheet 3.2g/ccm
Bv=Bv/1000;    %change units of fields from mT to T
Bs=Bs/1000;
yp='yes';
Ms   = 76; %70        %saturation magnetization
%Bc   = 0;        %circular field component [mT]

%% NEW
 rhy  = 56e-9;   %hydrodynamic radius [m]
 rco  = 12.5e-9;   %core radius [m]
 rho  = 5600;    %density [kg/m^3] from data sheet 3.2g/ccm
 Ms   = 70;      %saturation magnetization [emu/g]
 gam  = 1.3e9;   %gyromagnetic ratio [Hz/T]
 al   = 1;       %damping parameter
 K    = 4000;    %anisotropy constant [J/m^3]
%%

kT  = 1.38e-23*T;   %Boltzmann energy
Vhy = 4/3*pi*rhy^3; %NP hydrodynamic volume [m^3]
Vco = 4/3*pi*rco^3; %NP core volume [m^3]
mu  = rho*Ms*Vco   %magnetic moment [J/T]  would need rco=42e-9 to get the Micromod magnetization
%mu=Ms*10^(-18)*(rco/15e-9)^3;  %See uMod spec sheet
%mu=Ms*10^(-18); %See uMod spec sheet
%mu=3.2070e-18;
gm  = 6*visc*Vhy;   %drag coefficient [Nsm]
tE  = gm/2/kT;      %Einstein relaxation time [s]
YET = tE/sqrt(1+.21*(mu*Bv/kT)^2);

Dv=kT/gm; Du=mu/gm; %variables

%m  = initrand(0,1,N); %init mags randomly
%m   = zeros(N,3); m(:,3) = 1; % init mags in x
m = rand(N,3)*2-1; m = m./repmat(sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2),[1,3]); %size(m)
b = zeros(N,3); b(:,3)=Bv;  % make matrix b in z dir.
%calculate timesteps
if f==0; dt=10^-8; tf=5*tE;
else per=1/f; dt=per/tPts; tf=cycs*per;
end
t=0:dt:tf; disp(dt); AField=zeros(1,length(t));
M   = zeros(length(t),3);

%% Stochastic differential equation loop
in=1;

for t=0:dt:tf
    % save applied field
    AField(in)=Bv*cos(2*pi*f*t);
    h=randn(N,3); %stochastic term for fluctuations

    % Compute stochastic increments
    dW = sqrt(dt); %* randn(N,3);

    % Compute RK4 increments
    k1 = p(t, m, Du, b, Dv, Bs, f) * dt + q(m, Dv, h) * dW;
    k2 = p(t + 0.5 * dt, m + 0.5 * k1, Du, b, Dv, Bs, f) * dt + q(m + 0.5 * k1, Dv, h) * dW;
    k3 = p(t + 0.5 * dt, m + 0.5 * k2, Du, b, Dv, Bs, f) * dt + q(m + 0.5 * k2, Dv, h) * dW;
    k4 = p(t + dt, m + k3, Du, b, Dv, Bs, f) * dt + q(m + k3, Dv, h) * dW;

    % Update m using RK4 formula
    dm = (1/6) * (k1 + 2 * k2 + 2 * k3 + k4);
    m = m + dm;

    %normalize the magnitude
    nm = sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2);
    m(:,1)=m(:,1)./nm; m(:,2)=m(:,2)./nm; m(:,3)=m(:,3)./nm;

    % Save mean
    M(in,:)=sum(m);
    in=in+1;
end
t=0:dt:tf; %Return the time points

%DJ Plot the time-derivative of the magnetizaiton
%  load('my_colormap.mat'); colormap(my_colors);
%  figure; plot(t,M(:,1),'LineWidth',1,'Color',my_colors(1,:));
% hold on; plot(t,M(:,2),'LineWidth',1,'Color',my_colors(2,:));
%  hold on; plot(t,M(:,3),'LineWidth',1,'Color',my_colors(4,:));
%  legend('M_x','M_y (DC) ','M_z (AC)','Location','Southeast');
% %title('Magnetization versus time');
% title(['Viscosity: ', num2str(visc),' Pa-s; Time points: ', num2str(tPts), '; No. of NPs: ',num2str(N)]);
% xlabel('Time [s]'); ylabel('Magnetization');
%  set(gca,'FontWeight','Bold');
%
dt=t(2)-t(1);
tt=(t(1:end-1)+t(2:end))/2;
dMxdt=(M(2:end,1)-M(1:end-1,1))/dt;
dMydt=(M(2:end,2)-M(1:end-1,2))/dt;
dMzdt=(M(2:end,3)-M(1:end-1,3))/dt;
%  figure; plot(tt,dMxdt,'LineWidth',1,'Color',my_colors(1,:));
% hold on; plot(tt,dMydt,'LineWidth',1,'Color',my_colors(2,:));
% hold on; plot(tt,dMzdt,'LineWidth',1,'Color',my_colors(4,:));
% legend('dM_x/dt','dM_y/dt (DC) ','dM_z/dt (AC)','Location','Southeast');
% title('Derivative of magnetization versus time'); 
% title(['Viscosity: ', num2str(visc),' Pa-s; Time points: ', num2str(tPts), '; No. of NPs: ',num2str(N)]);
% xlabel('Time [s]'); ylabel('dM/dt'); %ylim([-4 4]*10^8);
% set(gca,'FontWeight','Bold');
dMdt(:,1)=dMxdt;
dMdt(:,2)=dMydt;
dMdt(:,3)=dMzdt;

end

% "helper" functions
function val=p(t, m, Du, b, Dv, Bs, f)
B=b*cos(2*pi*f*t); %drive field over time
B(:,1)=B(:,1)+Bs(1); B(:,2)=B(:,2)+Bs(2); B(:,3)=B(:,3)+Bs(3);
TB = Du*cross(cross(m,B),m); %mag torque
val = TB-2*m*Dv;
end

function val=q(m, Dv, h)
Ts=sqrt(2*Dv)*[h(:,2).*m(:,3)-h(:,3).*m(:,2),h(:,3).*m(:,1)-h(:,1).*m(:,3),h(:,1).*m(:,2)-h(:,2).*m(:,1)];
val = Ts;
end