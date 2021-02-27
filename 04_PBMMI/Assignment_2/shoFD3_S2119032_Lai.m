%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:Simple harmonic oscillator_lossy
%%%         Author:Nuolin Lai
%%%         Create Date:04/02/2021
%%%         Last modify date:06/02/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% parameters

SR = 44100;     % sample rate (Hz)
Tf = 1;         % duration (s)
f0 = 500;      % frequency (Hz)
u0 = 1;         % initial displacement
v0 = 0;         % initial velocity
T60 = 1;        % T60 values

k = 1/SR;               % time step
w0 = 2*pi*f0;           % angular frequency (rad./s)
Nf = floor(Tf*SR);      % total number of time steps

% lossy
sigma = 6*log(10)/T60; 
b = (2 - w0^2*k^2)/(1 + sigma*k);
a = (1 - sigma*k)/(1 + sigma*k);

% initialize

u2 = u0;                % set initial displacement
u1 = u0+k*v0;           % set second displacement

out = zeros(Nf,1);      % initialise output vector

% main loop

tic
for n=1:Nf
    u = b*u1 - a*u2;    % scheme update
    out(n) = u2;        % write output
    u2 = u1;            % shift state
    u1 = u;             % shift state
end
toc

% plot
M = length(out);
tax = [0:Nf-1]'*k;
subplot(2,1,1)
plot(tax,out,'k');
xlabel('t');
ylabel('u');
title('Simple Harmonic Oscillator-exact');
subplot(2,1,2);
spe = 20*log10(abs(fft(out)));
freq = (0:M-1)/M*SR;
semilogx(freq,spe);
xlabel('frequency (Hz)');
ylabel('Amplitude(dB)');
xlim([20 SR/2]);
% play sound

soundsc(out,SR)