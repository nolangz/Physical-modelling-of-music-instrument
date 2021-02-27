%--------------------------------------------
% PBMMI TUTORIAL 1 - part 2
% IIR FILTERS
% Alistair Carson 2021
% -------------------------------------------

clc
clear all
close all

% sample rate
Fs = 44100;

% filter parameters
f0 = 15e3;
Omega0 = 2*pi*f0/Fs;
r = 0.9;

% IR
N = 1024;                   % length
n = (0:N-1)';               % n-vector
h = 2*r.^n.*cos(Omega0*n);  % h[n]

% H(omega)
Hw = fft(h);
plot(2*(0:N-1)'/N, abs(Hw));
xlabel('Normalised frequency, \omega (\pi rad/s)'); ylabel('H(\omega)');
xlim([0 1]); grid on
title('FFT')

% filter coefficients
b = [1, -r*cos(Omega0)];
a = [1, -2*r*cos(Omega0), r^2];

% z-domain set-up
zmax = 2;
Nz = 100;
zv = linspace(-zmax, zmax, Nz);
[zr, zi] = meshgrid(zv, zv);
z = zr + 1i*zi;

% truncated FIR
Hz = zeros(size(z));
for n = 0:N-1
    Hz = Hz + h(n+1)*z.^(-n);
end

% analytic z-transform
Hz2 = (z./(z - r*exp(1i*Omega0)) + z./(z - r*exp(-1i*Omega0)));

% plot
plotZtransform(zr, zi, Hz, Hw, r);
title('Truncated')
plotZtransform(zr, zi, Hz2, Hw, r);
title('Analytic')

% Plots z transform, |z| = r, and FFT on same plot
function plotZtransform(zr, zi, Hz, Hw, r)
figure
surf(zr, zi, abs(Hz), 'FaceAlpha', 0.1);
xlabel('$\Re(z)$', 'Interpreter', 'latex')
ylabel('$\Im(z)$', 'Interpreter', 'latex')
zlabel('$|H(z)|$', 'Interpreter', 'latex')
zlim([0 12])
rotate3d on
hold on
Nfft = length(Hw);
theta = 2*pi*(0:Nfft-1)/Nfft;
plot(r*cos(theta), r*sin(theta), 'k', 'LineWidth', 2.0, 'LineStyle', '--');
plot3(cos(theta), sin(theta),abs(Hw), 'b', 'LineWidth', 2.0) 
legend('|H(z)|', '|z| = r', 'H(\omega)')
end