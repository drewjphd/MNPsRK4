%%%% Creates random vectors of [x,y,z] and length N with nonrandom part %%%%
%%%% DBR Sep 2011 %%%%

%inputs:
% number of triples (or particles) []
% partial alignment vector (x,y,z) format []

%outputs:
% a vector of Nx3 randomized or partially aligned triples

function A = initrand(N,b)

 A = zeros(N,3); 
 
 for ii=1:N
     r = 2*rand(1,3)-1+b;
 a=r/norm(r);
 
 A(ii,:) = a;
 end
