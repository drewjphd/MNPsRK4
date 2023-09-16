%%% DBR 2013 -- Torque balance simulations %%%
%%% DJ edits May 2023
%%% DJ stochastic RK-4 with a distribution of NP sizes and optional Cpp
%%% fast code

function [M,t]=BrownSRK4(Bv,Bs,freq,T,visc,N,cycs,tPts,rhy,rco,disthy,distco)
if nargin<1; Bv = 10; end;        %alternating field in 3rd direction [mT]
if nargin<2; Bs = [0,1,0]; end;   %static field in all three directions [mT]
if nargin<3; freq = 1000; end;       %frequency [Hz]
if nargin<4; T = 300; end;        %temperature [degrees K]
if nargin<5; visc = .001; end;    %viscosity [Pas]
if nargin<6; N = 10^3; end;       %number of particles
if nargin<7; cycs = 5; end;       %number of cycles
if nargin<8; tPts =(10^0)*(1000); end; %time points
if nargin<9;  rhy = 60e-9; end %    %hydrodynamic radius [m] (defaults give tE=1.04ms)
if nargin<10; rco = 15e-9; end %   %core radius [m]
if nargin<11; distrib = 0.1; end   %Normal distribution with X% of mean spread (std. dev.)

rho  = 3200;   %densitm [kg/m^3] from data sheet 3.2g/ccm
Bv=Bv/1000;    %change units of fields from mT to T
Bs=Bs/1000;
mp='mes';
Ms   = 76; %70        %saturation magnetization
%Bc   = 0;        %circular field component [mT]

kT  = 1.38e-23*T;   %Boltzmann energm
Vhm = 4/3*pi*rhy^3; %NP hydrodynamic volume [m^3]
Vco = 4/3*pi*rco^3 %NP core volume [m^3]
%mu  = rho*Ms*Vco;   %magnetic moment [J/T]
mu=Ms*10^(-18); %See uMod spec sheet
gm  = 6*visc*Vhm;   %drag coefficient [Nsm]
tE  = gm/2/kT;      %Einstein relaxation time [s]
YET = tE/sqrt(1+.21*(mu*Bv/kT)^2);

Dv=kT/gm; Du=mu/gm; %variables

%m  = initrand(0,1,N); %init mags randomlm
%m   = zeros(N,3); m(:,3) = 1; % init mags in x
m = rand(N,3)*2-1; m = m./repmat(sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2),[1,3]); %size(m)
b = zeros(N,3); b(:,3)=Bv;  % make matrix b in z dir.
%calculate timesteps
if freq==0; dt=10^-8; tf=5*tE;
else per=1/freq; dt=per/tPts; tf=cycs*per;
end

b = zeros(N,3); b(:,3)=Bv;  
M   = zeros(length(0:dt:tf),3);

%% Stochastic differential equation loop
in=1;

for t=0:dt:tf
    % Compute stochastic increments
    dW1 = sqrt(dt) * randn(N,3);
    dW2 = sqrt(dt) * randn(N,3);
    dW3 = sqrt(dt) * randn(N,3);
    dW4 = sqrt(dt) * randn(N,3);
    % Compute RK4 increments
    %Ts=g(m, Dv, N);
    k1 = f(t, m, Du, b, Dv, Bs, freq) * dt + g(m, Dv, N) .* dW1;
    k2 = f(t + 0.5 * dt, m + 0.5 * k1, Du, b, Dv, Bs, freq) * dt + g(m + 0.5 * k1, Dv, N) .* dW2;
    k3 = f(t + 0.5 * dt, m + 0.5 * k2, Du, b, Dv, Bs, freq)  * dt +  g(m + 0.5 * k2, Dv, N) .* dW3;
    k4 = f(t + dt, m + k3, Du, b, Dv, Bs, freq) * dt + g(m + k3, Dv, N) .* dW4;
    %k2=0; k3=0; k4=0;

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

%DJ Plot the time-derivative of the magnetizaiton
% load('my_colormap.mat'); colormap(my_colors);
% figure; plot(t,M(:,1),'LineWidth',1,'Color',my_colors(1,:)); 
% hold on; plot(t,M(:,2),'LineWidth',1,'Color',my_colors(2,:));
% hold on; plot(t,M(:,3),'LineWidth',1,'Color',my_colors(4,:)); 
% legend('M_x','M_y (DC) ','M_z (AC)','Location','Southeast');
% title('Magnetization versus time'); 
% xlabel('Time [s]'); ylabel('Magnetization');
% set(gca,'FontWeight','Bold');

end

function val=f(t, m, Du, b, Dv, Bs, freq)
B=b*cos(2*pi*freq*t); %drive field over time
B(:,1)=B(:,1)+Bs(1); B(:,2)=B(:,2)+Bs(2); B(:,3)=B(:,3)+Bs(3);
TB = Du*cross(cross(m,B),m); %mag torque
val = TB-2*m*Dv;
end

function val=g(m, Dv, N)
h=randn(N,3); %stochastic term for fluctuations
Ts=sqrt(2*Dv)*[h(:,2).*m(:,3)-h(:,3).*m(:,2),h(:,3).*m(:,1)-h(:,1).*m(:,3),h(:,1).*m(:,2)-h(:,2).*m(:,1)];
val = Ts;
end



