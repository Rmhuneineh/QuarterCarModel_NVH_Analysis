%% Clear ALL
clear all
close all
clc

%% Define Data
Ms = .4*2060; % [kg] Sprung Mass is mass distribution on the front axle multiplied by the vehicle's mass
Mu = .1*Ms; % [kg] Unsprung mass assumed 10% of sprung mass
Ks = 69516; % [N/m] Spring Stiffness on the AXLE not one wheel
Kt = 10*Ks; % [N/m] Unsprung mass spring stiffness assumed 10 times KS
Cs = sqrt(Ks*Ms)/2; % [Ns/m] First attempt assumes the optimal value for Cs

%% State Space Representation
M = [Ms 0; 0 Mu];
K = [Ks -Ks; -Ks Ks+Kt];
Cc = [Cs -Cs; -Cs Cs];
In = [0; Kt];
A = [zeros(2, 2) eye(2);-M\K -M\Cc];
B = [zeros(2, 1); M\In];
To = [1 0];
C = [-To*(M\K) -To*(M\Cc)];
D = To*(M\In);

%% Bode Plots
G = ss(A, B, C, D);
omega = 2*pi*(0.01:0.1:80);
[G_mag, ~] = bode(G, omega);
[Ny, Nu, lw] = size(G_mag);
G_mag = reshape(G_mag, [Ny*Nu, lw]).';

num = [80.03 989 0.02108];
den = [1 78.92 2412 5614];
H_2631 = tf(num, den);
GH = G*H_2631;
[GH_mag, ~] = bode(GH, omega);
[Ny, Nu, lw] = size(GH_mag);
GH_mag = reshape(GH_mag, [Ny*Nu, lw]).';

figure(1)
loglog(omega/(2*pi), G_mag)
hold on
loglog(omega/(2*pi), GH_mag);
grid on

%% Power Spectral Density
c = 6.4*1e-7; % [m^2.cycles/m] - Road Grade B
V = 70/3.6; % [m/s] - Vehicle Velocity
f = (0.1:0.1:80)'; % [Hz]
Sin = c*V*f.^(-2);
Sout = G_mag.^2.*Sin;
Socc = GH_mag.^2.*Sin;

figure(2)
loglog(f, Sin);
hold on
loglog(f, Sout);
loglog(f, Socc);
grid on

%% Optimizing Cs Relative To Mean Square Value
Cs = 1e2:1:1e4;
MSV_qs = zeros(1, length(Cs));
MSV_qocc = zeros(1, length(Cs));
for i = 1:length(Cs)
    Cc = [Cs(i) -Cs(i); -Cs(i) Cs(i)];
    A = [zeros(2, 2) eye(2);-inv(M)*K -inv(M)*Cc];
    C = [-To*inv(M)*K -To*inv(M)*Cc];
    G = ss(A, B, C, D);
    omega = 2*pi*(0.1:0.1:80);
    [G_mag, ~] = bode(G, omega);
    [Ny, Nu, lw] = size(G_mag);
    G_mag = reshape(G_mag, [Ny*Nu, lw]).';
    GH = G*H_2631;
    [GH_mag, ~] = bode(GH, omega);
    [Ny, Nu, lw] = size(GH_mag);
    GH_mag = reshape(GH_mag, [Ny*Nu, lw]).';
    Sout = G_mag.^2.*Sin;
    Socc = GH_mag.^2.*Sin;
    MSV_qs(i) = trapz(f, Sout);
    MSV_qocc(i) = trapz(f, Socc);
end

figure(3)
loglog(Cs, MSV_qs)
hold on
loglog(Cs, MSV_qocc)
grid on