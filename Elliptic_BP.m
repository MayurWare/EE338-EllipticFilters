%Elliptic Lowpass filter transfer function and plots
[z, p, H0, B, A] = ellipap2(4, 1.411, 16.48);
syms sL Omega_L;
H_LPF = (0.8501*(1 + 0.8834*sL^2)*(1+0.3195*sL^2))/((1+0.0620*sL+sL^2)*(1+1.1316*sL+1.5999*sL^2));
H_LPF_f = subs(H_LPF, sL, 1i * Omega_L);
[ns, ds] = numden(H_LPF);
nsl = sym2poly(ns);
dsl = sym2poly(ds);
k = ds(1);
nsl = nsl / k;
dsl = dsl / k;
disp(nsl);
disp(dsl);
figure();
fplot(abs(H_LPF_f));

%Analog lowpass to bandpass frequency transformation
syms s Omega;
Omega0 = 0.8265;
B0 = 0.449;
H_BPF = subs(H_LPF, sL, (s^2 + Omega0^2) / (B0 * s));
[ns, ds] = numden(H_BPF);
ns = sym2poly(ns);
ds = sym2poly(ds);
k = ds(1);
ns = ns / k;
ds = ds / k;
disp(ns);
disp(ds);
H_BPF_f = subs(H_BPF, s, 1i * Omega);
figure();
fplot(abs(H_BPF_f));

%Analog to z bilinear transformation
syms z;
Hz = subs(H_BPF, s, (z-1)/(z+1));
[nz, dz] = numden(Hz);
nz = sym2poly(nz);
dz = sym2poly(dz);
k = dz(1);
nz = nz / k;
dz = dz / k;
disp(nz);
disp(dz);
[H,f] = freqz(nz,dz,1024*1024, 540e3);
plot(f,abs(H));
fvtool(nz, dz, 'Analysis', 'Phase');

z = roots(nz);
p = roots(dz);
zplane(z,p);
grid
