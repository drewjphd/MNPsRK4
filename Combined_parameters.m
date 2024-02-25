function [t0,tN,tB,xi0,sig,A,B1,B2]=Combined_parameters(Bv,f)

if nargin==0; Bv=0; f=0; end
    
T    = 280; %293;    %temperature [K]
visc = .001;    %viscosity [Pas]
rhy  = 56e-9;   %hydrodynamic radius [m]
rco  = 12.5e-9;   %core radius [m]
rho  = 5600;    %density [kg/m^3] from data sheet 3.2g/ccm
Ms   = 70;      %saturation magnetization [emu/g]
gam  = 1.3e9;   %gyromagnetic ratio [Hz/T]
al   = 1;       %damping parameter
K    = 4000;    %anisotropy constant [J/m^3]

kT   = 1.38e-23*T;   %Boltzmann energy  
Vh   = 4/3*pi*rhy^3; %NP hydrodynamic volume [m^3]
Vc   = 4/3*pi*rco^3; %NP core volume [m^3]
mu   = rho*Ms*Vc;    %magnetic moment [J/T]
xi0  = mu*Bv/kT;     %unitless field []
tB   = 3*visc*Vh/kT; %Einstein relaxation time [s]
sig  = K*Vc/kT;      %unitless anisotropy

t0   = mu/(2*gam*kT)*(1+al^2)/al; %measurement time [s]

  if sig<1 %limits from Garcia Palacios
  tN=t0*(1-2/5*sig+48/875*sig^2)^(-1);
  else
  tN=t0*exp(sig)/2*sqrt(pi./sig);
  end
  
A=xi0/f/tB;

B1=xi0/f/t0;
B2=sig/f/t0;
