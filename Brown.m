%%% YP 2016 -- Torque balance simulations %%%
%compensating the timepoints effect 7/22/2016/yipeng
%%% 1/13/2017 use varying timesteps for simulation to avoid distortion
function [R,Ha,Hb]=Brown(Bv,T,f,visc,mu)
%Ha Hb is H2 H4 or H3 H5 (Bs=0)
bb=Bv/1000;
yp='yes'; %'yes' means plotting graphs
Bs   = 0;         %static field [mT] %if Bs is zero, do inline simulation
cycs = 5;      %number of cycles
%f    = 10;    %frequency [Hz]
%Bv   = 10;     %oscillating field amp [mT]
% tPts = 2^10;    %time points
%T    = 293;    %temperature [K]
rhy  = 57e-9;  %hydrodynamic radius [m]
rco  = 25e-9;  %core radius [m]
rho  = 3200;   %density [kg/m^3] from data sheet 3.2g/ccm

N    = 10^4;      %number of particles
% visc = .001;   %viscosity [Pas]

%Bc   = 0;        %circular field component [mT]
kT  = 1.38e-23*T;   %Boltzmann energy
Vhy = 4/3*pi*rhy^3; %NP hydrodynamic volume [m^3]
% Vco = 4/3*pi*rco^3; %NP core volume [m^3]
% Ms   = 70;        %saturation magnetization [emu/g]
% mu  = rho*Ms*Vco;   %magnetic moment [J/T]   %%%8-31-15 this way of
%   calculating mu seems to be inaccurate since Ms is not accurate
%mu=7.6*10^(-17); %data from Micromod particle data sheet
gm  = 6*visc*Vhy;   %drag coefficient [Nsm]
tE  = gm/2/kT;      %Einstein relaxation time [s]
%YET = tE/sqrt(1+.21*(mu*bb/kT)^2);

Dv=kT/gm; Du=mu/gm; %variables
tPts=max(floor(Du*Bv/f/800),500); %change timepoints based on 0.2% error
%tPts=10*tPts
tPts=10^3;
%m  = initrand(0,1,N); %init mags randomly
m   = zeros(N,3); m(:,1) = 1; %init mags in x
b   = zeros(N,3); b(:,3)=Bv/1000; %make matrix b in z dir.
%calculate timesteps
if f==0; dt=10^-5; tf=5*tE;
else per=1/f; dt=per/tPts; tf=cycs*per;
end
t=0:dt:tf;
M   = zeros(length(t),3);

%% Stochastic differential equation loop
in=1;
for i=0:dt:tf
    h=randn(N,3); %stochastic term for fluctuations
    B=b*cos(2*pi*f*i); %drive field over time
    B(:,1)=Bs/1000;
    
    %alpha=mu*Bv*cos(2*pi*f*i)/1000/kT; ML=coth(alpha)-1/alpha;
    %sm=.2;%std(m(:,3));
    %good=find(m(:,3)<ML+sm & m(:,3)>ML-sm) ;
    
    % mm(in,:,:)=m; %for single particle trajectories
    %M(in,:)=mean(m(good,:));
    M(in,:)=mean(m);
    if in==2600
        subplot(1,3,1)
        plot(M(:,3));
        subplot(1,3,2);
        plot(M(:,1));
        subplot(1,3,3);
        plot(M(:,2));
    end
    TB=Du*cross(cross(m,B),m); %mag torque
    Ts=sqrt(2*Dv)*cross(h,m); %stochastic torque
    
    %stochastic diffeq Stratanovich with drift term
    dm = (TB-2*m*Dv)*dt + Ts*sqrt(dt);
    m = m + dm;
    
    %normalize the magnitude
    nm = sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2);
    m(:,1)=m(:,1)./nm; m(:,2)=m(:,2)./nm; m(:,3)=m(:,3)./nm;
    
    %good(:,in)=find(m(:,3) < m(:,3)+2*std(m));
    
    
    in=in+1;
end;

%% Plotting and harmonic analysis
if f~=0
    lastper=2*tPts+1:3*tPts; H=zeros(10,3);
    pM=(fft(M(lastper,:)))/tPts; %compensating the time points effect
    for i =1:10; H(i,:)=pM(i+1,:)*i; end %find the derivative harmonics
    
    if Bs~=0
        Ha=H(2,1); Hb=H(4,1); %second and fourth harmonics
    else
        Ha=H(3,3); Hb=H(5,3); %third and fifth harmonics
    end
    R=abs(Hb./Ha); %MSB signal
    
    fMz=fft(M(lastper,3));
    
    %R42=H(4,1)/H(2,1); %sMSB signal
    if strcmp(yp,'yes');
        figure;
        %    figuresize(8,4,'inches' )
        subplot(1,7,1:4); plot(t*1000,M)
        xlabel('Time (ms)'); ylabel('Normalized mean magnetization')
        legend('Mx','My','Mz','Location','NorthEast');
        legend('boxoff')
        subplot(1,7,6:7); bar(abs(H))
        xlabel('Harmonic #');  ylabel('Signal Magnitude');
        legend('Hx','Hy','Hz','Location','NorthEast');
        legend('boxoff')
        xlim([0 10])
    end
    MM=max(M(:,3)); alpha=mu*bb/kT; ML=coth(alpha)-1/alpha;
else
    MM=M(length(M),3);
    alpha=mu*bb/kT; ML=coth(alpha)-1/alpha;
    H3=0; H5=0; R53=0; %MSB signall
    
end

