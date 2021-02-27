%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:BTB 2# raw code
%%%         Author:Nuolin Lai
%%%         Create Date:29/01/2020
%%%         Last modify date:30/01/2020
%%%         1.decay_stretching_tuning
%%%         2.rests_and_end_notes
%%%         3.guitar body impulse response convolution
%%%         4.moving pick
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

clc;
clear all;
close all;

%set governing paramenter

Fs = 44.1e3;                            %sampling rate in Hz
T60 = 9;                                %duration of simulation in seconds
tau = T60/log(1000);                    % calculate tau
f0 =73;                                 %fundamental frequency of the string in Hz


rho =  exp(-1/(f0*tau))/abs(cos(pi*f0/Fs)); %loss parameter rho
R = 0.95;                               % dynamics parameter R
S = 0.5;                                % the decay stretching factor
P1 = Fs/f0;                             % the real value for the period of the first partial
Pa = S;                                 %phase delay
%omega0 = 2*pi*f0;
% if f0 <= Fs/40
%     Pa = S;
% else
%     Pa = 1/(tan(-S*sin(omega0/Fs)/((1-S)-S*cos(omega0/Fs)))*omega0/Fs);
% end
Nexact = P1-Pa;                         %ideal number of samples
M = T60*Fs;                             %the length of output 
N = floor(Nexact);                      %number of sample in delay line
P = Nexact-N;                           %fractional delay
C = (1-P)/(1+P);                        %Gain coefficient in Hc

% input and ouput vector initialise
v = -1 + (1+1)*rand(N+1,1);             %create white noise from -1 to 1
y = zeros(M,1);                         %output of the algorithm

%%
%moving pick based on comb filter on white noise
Mu= 1/2;                                %position of pick
delay_s = floor(Mu*N);                  %delay sample of comb filter
dlinebuf = zeros(delay_s,1);            %delay line buffer

%circular buffer to create comb filter
for n = 1:N+1
    vp(n)   = -dlinebuf(mod(n-1,delay_s)+1)+v(n);
    %update buffer
    dlinebuf(mod((n+delay_s-2),delay_s)+1) = v(n);
end

%%
%implement dynamics filter

x1 = 0; %initialise x1

%read from white noise
for n = 0:N
    x0 = (1-R)*vp(n+1);
    y(n+1) = x0+R*x1;
    x1 = x0;
end

%decay stretching and tuning
yp1 = 0;
for n = N+1:M-1
    yp0= y(n-N+1)*C+y(n-N)+yp1*(-C) ;    %output of Hc
    y(n+1) = rho*((1-S)*yp0+S*yp1);      %output of Ha
    yp1 = yp0;
end

%%
%convolution
IR_guitar = audioread('Taylor814.wav');
y = conv(y,IR_guitar);
L = length(y);

%%
%Normalise and plot
y = y/max(y);
t = (0:L-1)/Fs;
subplot(2,1,1);
plot(t,y);
xlabel('Time (s)');
ylabel('Amplitude');
title('Waveform')
subplot(2,1,2);
spe = 20*log10(abs(fft(y)));
freq = (0:L-1)/L*Fs;
semilogx(freq,spe);
xlabel('frequency (Hz)');
ylabel('Amplitude(dB)');
title('Spectrum')
ylim([-20 100]);
xlim([20 Fs/2]);
xline(f0,'--')

sound(y,Fs);