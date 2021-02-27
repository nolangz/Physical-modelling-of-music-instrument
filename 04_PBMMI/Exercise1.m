%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:Exercise1
%%%         Author:Nuolin Lai
%%%         Create Date:18/01/2020
%%%         Last modify date:18/01/2020
%%%         FIR filter 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;

%import audio
[x,Fs]   = audioread('birchcanoe.wav');
N =length(x);

%change x to N-length Kronecker delta function 
x = [1;zeros(N-1,1)]

%create impulse response vector h
Nb = 2;
n_coef = Nb+1;
h = [1,1,1];

%frequency vector
fx = (0:N-1)';
fy = (0:N-1+2*Nb)';

%create output signal
y = h(1)*[x;0;0]+h(2)*[0;x;0]+h(3)*[0;0;x];

%DFT
XF = fft(x);
YF = fft(y);
L1 = length(XF);
plot(fy,20*log10(abs(YF)));

T = fvtool(h,1);

%play the sound
sound(y,Fs)
