%Bandpass Elliptical Filter Functions
function [z,p,H0,B,A] = ellipap2(N,Ap,As)
if nargin==0, help ellipap2; return; end
%Defining Passband Gain
Gp = 10^(-Ap/20); 
%Ripple Factors
ep = sqrt(10^(Ap/10) - 1);
es = sqrt(10^(As/10) - 1);
k1 = ep/es;
k = 0.8723;
%L is the number of second-order sections
L = floor(N/2); r = mod(N,2); 
%Zeros of elliptic rational function
i = (1:L)'; ui = (2*i-1)/N; zeta_i = cde(ui,k); 
%Solving for Poles and Zeros
z = j./(k*zeta_i); 
v0 = -j*asne(j/ep, k1)/N; 
p = j*cde(ui-j*v0, k); 
p0 = j*sne(j*v0, k); 
B = [ones(L,1), -2*real(1./z), abs(1./z).^2]; 
A = [ones(L,1), -2*real(1./p), abs(1./p).^2];
if r==0, 
B = [Gp, 0, 0; B]; A = [1, 0, 0; A];
else
B = [1, 0, 0; B]; A = [1, -real(1/p0), 0; A];
end
%Conjugate Zeros
z = cplxpair([z; conj(z)]); 
%Conjugate Poles
p = cplxpair([p; conj(p)]);
if r==1, p = [p; p0]; end 
H0 = Gp^(1-r); 
end