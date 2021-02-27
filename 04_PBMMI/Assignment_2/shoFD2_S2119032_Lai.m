%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:Simple harmonic oscillator(exact)
%%%         Author:Nuolin Lai
%%%         Create Date:04/02/2021
%%%         Last modify date:05/02/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% parameters

SR = 44100;     % sample rate (Hz)
Tf = 1;         % duration (s)
f0 = 15e3;     % frequency (Hz)
u0 = 1;         % initial displacement
v0 = 0;         % initial velocity


k = 1/SR;               % time step
w0 = 2*pi*f0;           % angular frequency (rad./s)
Nf = floor(Tf*SR);      % total number of time steps

alpha = 2/w0^2/k^2-cos(w0*k)/(1-cos(w0*k));  % the scheme free parameter

b = (2-alpha*w0^2*k^2)/(1 + ((1-alpha)*w0^2*k^2)/2); %coefficient of u1

be = 2*cos(w0*k); %coefficient of u1(exact)

%assert
if(k > 2/(w0*sqrt(2*alpha-1)) && alpha >= 0.5)
    error('Stability condition violated')
end

% initialize
u2 = u0;                        % set initial displacement
u1 = u0+k*v0+k^2/2*(-w0^2*u0);  % set second displacement %3rd order

u2e = u0;                       % set initial displacement(exact)
u1e = u0+k*v0+k^2/2*(-w0^2*u0); % set second displacement(exact) %3rd order

%initialise
out = zeros(Nf,1);      % output vector
outex = zeros(Nf,1);    % output vector
er = zeros(Nf,1);       % error vector

% main loop
tic
for n=1:Nf
    u = b*u1 - u2;           % scheme update
    out(n) = u2;             % write output
    u2 = u1;                 % shift state
    u1 = u;                  % shift state
    ue = be*u1e - u2e;       % scheme update
    outex(n) = u2e;          % write output
    u2e = u1e;               % shift state
    u1e = ue;                % shift state
    er(n) = out(n)-outex(n); % calculate error
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
speex = 20*log10(abs(fft(outex)));
freq = (0:M-1)/M*SR;
semilogx(freq,spe,freq,speex);
xlabel('frequency (Hz)');
ylabel('Amplitude(dB)');
xlim([20 SR/2]);
% play sound

soundsc(out,SR)
