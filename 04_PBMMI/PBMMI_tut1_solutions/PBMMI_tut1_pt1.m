%------------------------------------------------
% PBMMI TUTORIAL 1
% Part 1: FIR filters and convolution
% ----------------------------------------------
clc
clear all
close all

% input
load handel; x=y; clear y;
Nx = length(x);

% kronecker delta
%x = [1; zeros(Nx-1, 1)];

% IR (window function)
Nb = 100;                                % filter order
h = ones(Nb+1, 1);                       % rectangular
%h = 0.5 - 0.5*cos(2*pi*(0:Nb)'/(Nb+1)); % hann

% Filter and compare timings
tic
for n = 1:10
    y = myconv(h, x);
end
toc
tic
for n = 1:10
    y2 = conv(h, x);
end
toc
tic
for n = 1:10
    y3 = myfastconv(h, x);
end
toc

% Plotting
plot_spectra(x, y, Fs);


function plot_spectra(x, y, Fs)
Nx = length(x);
Ny = length(y);
fx = Fs*(0:Nx-1)'/Nx; 
fy = Fs*(0:Ny-1)'/Ny;
XF = fft(x); 
YF = fft(y);
plot(fx,20*log10(abs(XF)),fy,20*log10(abs(YF))); 
xlim([0, Fs/2]);
xlabel('Frequency (Hz)'); ylabel('dB');
end