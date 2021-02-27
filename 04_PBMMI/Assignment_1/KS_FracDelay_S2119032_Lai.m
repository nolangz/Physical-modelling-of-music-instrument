%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:KS_FracDelay
%%%         Author:Nuolin Lai
%%%         Create Date:25/01/2020
%%%         Last modify date:28/01/2020
%%%         Part 2 of the assignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;

%set governing paramenter

Fs = 44.1e3; %sampling rate in Hz
dur = 5; %duration of simulation in seconds
f0 = 440; %fundamental frequency of the string in Hz
rho =  0.99; %loss parameter rho
R = 0.95; % dynamics parameter R
Nexact = Fs/f0-0.5;%ideal number of samples

%Calculate the derived parameter

M = dur*Fs; % the length of output 
N = floor(Nexact); %number of sample in delay line
P = Nexact-N;%fractional delay
C = (1-P)/(1+P);
% vector initialise
rng(0);
v = -1 + (1+1)*rand(N+1,1); %create white noise from -1 to 1
rng(0);
y = zeros(M,1); %output of the algorithm

%implement dynamics filter
x1 = 0; %initialise x1

%read from white noise
for n = 0:N
    x0 = (1-R)*v(n+1);
    y(n+1) = x0+R*x1;
    x1 = x0;
end

%KS main implement 1
% tic
% yp(N+1) = 0;
% for n = N+1:M-1
%     yp(n+1) = y(n+1-N)*C+y(n-N)+yp(n)*(-C) ;
%     y(n+1) = rho/2*(yp(n+1)+yp(n));
% end
% t1=toc;

%KS main implement 2
tic
yp1 = 0;
for n = N+1:M-1
    yp0= y(n-N+1)*C+y(n-N)+yp1*(-C) ;
    y(n+1) = rho/2*(yp0+yp1);
    yp1 = yp0;
end
t2=toc;

% t1 = 0.0221
% t2 = 0.0033

t = (0:M-1)/Fs;
subplot(2,1,1);
plot(t,y);
%xlabel
xlabel('Time (s)');
ylabel('Amplitude');
title('Waveform')
subplot(2,1,2);
spe = 20*log10(abs(fft(y)));
freq = (0:M-1)/M*Fs;
semilogx(freq,spe);
%xlabel
xlabel('frequency (Hz)');
ylabel('Amplitude(dB)');
title('Spectrum')
ylim([-20 100]);
xlim([20 Fs/2]);
xline(f0,'--')

sound(y,Fs)