%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:KS_Basic
%%%         Author:Nuolin Lai
%%%         Create Date:18/01/2021
%%%         Last modify date:18/01/2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;

%% set governing paramenter

Fs = 44.1e3; %sampling rate in Hz
dur = 5; %duration of simulation in seconds
f0 = 73; %fundamental frequency of the string in Hz
rho =  0.998; %loss parameter rho
R = 0.95; % dynamics parameter R

%% Calculate the derived parameter

M = dur*Fs; % the length of output 
N = ceil(Fs/f0-1/2); %delay line length

%%  vector initialise
rng(0);
v = -1 + (1+1)*rand(N+1,1); %create white noise from -1 to 1
rng(0);
y = zeros(M,1); %output of the algorithm

%% implement dynamics filter
x1 = 0; %initialise x1

%read from white noise
for n = 0:N
    x0 = (1-R)*v(n+1);
    y(n+1) = x0+R*x1;
    x1 = x0;
end

%implement main KS algorithm
for n = N+1:M-1
    y(n+1) = rho/2*(y(n-N+1)+y(n-N));
end
%% plot waveform and spectrum
t = (0:M-1)/Fs;
subplot(2,1,1);
plot(t,y);
xlabel('Time (s)');
ylabel('Amplitude');
title('Waveform')
subplot(2,1,2);
spe = 20*log10(abs(fft(y)));
freq = (0:M-1)/M*Fs;
semilogx(freq,spe);
xlabel('frequency (Hz)');
ylabel('Amplitude(dB)');
title('Spectrum')
ylim([-20 100]);
xlim([20 Fs/2]);
xline(f0,'--');

%% play
sound(y,Fs)