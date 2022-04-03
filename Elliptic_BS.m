%Elliptic Lowpass filter transfer function and plots
[z, p, H0, B, A] = ellipap2(3, 1.411, 16.48);
syms sL Omega_L;
H_LPF = (1 + 0.6294*sL^2)/((1+0.2306*sL+0.9994*sL^2)*(1+1.6046*sL));
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
Omega0 = 0.7496;
B0 = 0.389;
H_BSF = subs(H_LPF, sL, (B0 * s)/(s^2 + Omega0^2));
[ns, ds] = numden(H_BSF);
ns = sym2poly(ns);
ds = sym2poly(ds);
k = ds(1);
ns = ns / k;
ds = ds / k;
disp(ns);
disp(ds);
H_BPF_f = subs(H_BSF, s, 1i * Omega);
figure();
fplot(abs(H_BPF_f));

%Analog to z bilinear transformation
syms z;
Hz = subs(H_BSF, s, (z-1)/(z+1));
[nz, dz] = numden(Hz);
nz = sym2poly(nz);
dz = sym2poly(dz);
k = dz(1);
nz = nz / k;
dz = dz / k;
disp(nz);
disp(dz);
[H,f] = freqz(nz,dz,1024*1024, 400e3);
plot(f,abs(H));
fvtool(nz, dz, 'Analysis', 'Phase');

z = roots(nz);
p = roots(dz);
zplane(z,p);
grid
