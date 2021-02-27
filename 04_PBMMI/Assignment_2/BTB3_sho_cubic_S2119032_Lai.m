%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:Simple harmonic oscillator_non-linear
%%%         Author:Nuolin Lai
%%%         Create Date:04/02/2021
%%%         Last modify date:07/02/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% parameters

SR = 44100;     % sample rate (Hz)
Tf = 3;         % duration (s)
f0 = 500;       % frequency (Hz)
u0 = 1.185;     % initial displacement
v0 = 0;         % initial velocity

alpha = 0.5;    % the scheme free parameter

k = 1/SR;               % time step
       
w1 = sqrt(2*f0*pi);     % angular frequency (rad./s)
Nf = floor(Tf*SR);      % total number of time steps


%assert
if(k > 2/(w1*sqrt(2*alpha-1)) && alpha >= 0.5)
    error('Stability condition violated')
end

b = (2-alpha*w1^2*k^2)/(1 + ((1-alpha)*w1^2*k^2)/2); %coefficient of u1

% initialize

u2 = u0;                                   % set initial displacement
u1 = u0+k*v0+k^2/2*(-w1^4*u0^3);           % set second displacement(3rd order)
out = zeros(Nf,1);                         % output vector

% main loop

tic
for n=1:Nf
    u = 2/(1+k^2*w1^4*u1^2/2)*u1 - u2;        % scheme update
    out(n) = u2;                              % write output 
    u2 = u1;                                  % shift state
    u1 = u;                                   % shift state
end
toc

% plot
M = length(out);
tax = [0:Nf-1]'*k;
subplot(2,1,1)
plot(tax,out,'k');
xlabel('t');
ylabel('u');
title('Simple Harmonic Oscillator-nonlinear');
subplot(2,1,2);
spe = 20*log10(abs(fft(out)));
freq = (0:M-1)/M*SR;
semilogx(freq,spe);
xlabel('frequency (Hz)');
ylabel('Amplitude(dB)');
xlim([20 SR/2]);
% play sound

soundsc(out,SR)
