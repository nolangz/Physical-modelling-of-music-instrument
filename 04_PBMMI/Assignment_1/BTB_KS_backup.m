%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:function for single notes
%%%         Author:Nuolin Lai
%%%         Create Date:18/01/2021
%%%         Last modify date:18/01/2021
%%%         Input BPM,syncopated note type(1,1/2,1/4),F0,Fs,R,S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

%set global paramenter
Fs  = 44.1e3;         %Sampling rate in Hz
T60 = 5;              %Duration of simulation in seconds
r   = T60/log(1000);  %Calculate tau
f0  =75;              %Fundamental frequency of the string in Hz
f0s = f0/2;           %End frequency of slur
BPM = 150;            %BPM
Bt = 60/BPM;          %beat duration in seconds.
rho =  exp(-1/(f0*r))/abs(cos(pi*f0/Fs)); %loss parameter rho
R = 0.95;                                 % dynamics parameter R
S = 0.5;                                  % the decay stretching factor
P1 = Fs/f0;   % the real value for the period of the first partial(start)
Ps = Fs/f0s;  % the real value for the period of the first partial(end)
Pa = S;       % phase delay for low frequency
%omega0 = 2*pi*f0; 
% if f0 <= Fs/40
%     Pa = S;
% else
%     Pa = 1/(tan(-S*sin(omega0/Fs)/((1-S)-S*cos(omega0/Fs)))*omega0/Fs);
% end
Nexact = P1-Pa;         %ideal number of samples for original pitch
Nexacts = Ps-Pa;        %ideal number of samples for end pitch
N = floor(Nexact);      %number of sample in delay line(start)
N2 = floor(Nexacts);    %number of sample in delay line(end)

M = T60*Fs;             %the length of output 
Gt = 2*Bt;              %time for glissandis
NG = Gt*Fs;             %number of sample in glissandis process
Ns = N*ones((M-N-1),1); %vector of delay line 
Gs = linspace(N,N2,NG); %Delay sample during glissandis process
Ns(1:NG)=Gs;
Ns = floor(Ns);

fv = 5; %frequency of vibrato
Vs = min(N,N2)+abs(N-N2)/2+abs(N-N2)/2*sin(2*pi*fv*(0:NG-1)/Fs); % delay vector of vibrato

P = Nexact-N;%fractional delay
C = (1-P)/(1+P);
S = 0.5;
% vector initialise
v = -1 + (1+1)*rand(N+1,1); %create white noise from -1 to 1
y1 = zeros(M,1); %output of the algorithm
% 
%implement dynamics filter
x1 = 0; %initialise x1

%read from white noise
% for n = 0:N
%     x0 = (1-R)*v(n+1);
%     y(n+1) = x0+R*x1;
%     x1 = x0;
% end
% tic;
% %decay stretching and tuning
% yp1 = 0;
% 
% for n = N+1:M-1
%     yp0= y(n-N+1)*C+y(n-N)+yp1*(-C) ;
%     y(n+1) = rho*((1-S)*yp0+S*yp1);
%     yp1 = yp0;
% end
% 
x1 = 0; %initialise x1
%read from white noise
for n = 0:N
    x0 = (1-R)*v(n+1);
    y(n+1) = x0+R*x1;
    x1 = x0;
end

% Slur
tic
yp1 = 0;
for n = N+1:M/5
    yp0= y(n-N+1)*C+y(n-N)+yp1*(-C) ;
    y(n+1) = rho*((1-S)*yp0+S*yp1);
    yp1 = yp0;
end

for n = M/5+1:M-1
    yp0= y(n-N2+1)*C+y(n-N2)+yp1*(-C) ;
    y(n+1) = rho*((1-S)*yp0+S*yp1);
    yp1 = yp0;
end
%
% x1 = 0; %initialise x1
% %read from white noise
% for n = 0:N
%     x0 = (1-R)*v(n+1);
%     y2(n+1) = x0+R*x1;
%     x1 = x0;
% end
% 
% %glissandis
% yp1 = 0;
% for n = N+1:M-1
%     yp0= y2(n-Ns(n-N)+1)*C+y2(n-Ns(n-N))+yp1*(-C);
%     y2(n+1) = rho*((1-S)*yp0+S*yp1);
%     yp1 = yp0;
% end
% IR_guitar = audioread('Taylor814.wav');
% y2 = conv(y2,IR_guitar);
L = length(y);

%Normalise
y = y/max(y);
t = (0:L-1)/Fs;
subplot(2,1,1);
plot(t,y);
%xlabel
xlabel('Time (s)');
ylabel('Amplitude');
title('Waveform')
subplot(2,1,2);
spe = 20*log10(abs(fft(y)));
freq = (0:L-1)/L*Fs;
semilogx(freq,spe);
%xlabel
xlabel('frequency (Hz)');
ylabel('Amplitude(dB)');
title('Spectrum')
ylim([-20 100]);
xlim([20 Fs/2]);
xline(f0,'--')

%
sound(y,Fs)