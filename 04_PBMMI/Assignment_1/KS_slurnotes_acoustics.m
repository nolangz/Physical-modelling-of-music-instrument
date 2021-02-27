%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:KS_slurnotes_acoustics
%%%         Author:Nuolin Lai
%%%         Create Date:28/01/2021
%%%         Last modify date:18/01/202
%%%         Input BPM,syncopated note type(1,1/2,1/4),fundamental freqeuncy
%%%         Sampling frequency,Dynamics parameter,decay stretching factor
%%%         BPM = 60~120 (suggestion)
%%%         syn1,2 = 4,2,1,1/2,1/4,1/8....
%%%         f0,1   = frequency of notes
%%%         Fs: sampling rate
%%%         R（0，1）: dynamics parameter R
%%%         S（0，1）: the decay stretching factor
%%%         Mu (0,1): pick position (0,1) 1/2 is in the middle
%%%         example y= KS_slurnotes_acoustics(60,2,2,100,150,44100,0.95,0.9,0.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function y = KS_slurnotes_acoustics(BPM,syn1,syn2,f0,f1,Fs,R,S,Mu);
%% set global paramenter
Bt  = 60/BPM;          % beat duration in seconds.
factor = 2;            % correct factor between feeling duration and T60
T60 = Bt*(syn1+syn2)*factor;  % duration for sum of two syncopated note
M1  = floor(Bt*syn1*Fs);      % duration for the first syncopated note
Fs  = 44.1e3;          % sampling rate in Hz
tau = T60/log(1000);   % calculate tau         
rho =  exp(-1/(f0*tau))/abs(cos(pi*f0/Fs)); %loss parameter rho
P1  = Fs/f0;           % the real value for the period of the first partial
P12  = Fs/f1;          % the real value for the period of the first partial
Pa = S;
Nexact = P1-Pa;             % ideal number of samples
Nexact2 = P12-Pa;           % ideal number of samples
M = floor(T60*Fs);          % the length of output 
N = floor(Nexact); 
N2 = floor(Nexact2);        % number of sample in delay line
P = Nexact-N;               % fractional delay
C = (1-P)/(1+P);            % dynamics parameter R
S = 0.5;                    % the decay stretching factor
% vector initialise
v = -1 + (1+1)*rand(N+1,1); % create white noise from -1 to 1
y = zeros(M,1);             % output of the algorithm

%% moving pick based on comb filter on white noise
delay_s = floor(Mu*N);                  %delay sample of comb filter
dlinebuf = zeros(delay_s,1);            %delay line buffer

%circular buffer to create comb filter
for n = 1:N+1
    vp(n)   = -dlinebuf(mod(n-1,delay_s)+1)+v(n);
    %update buffer
    dlinebuf(mod((n+delay_s-2),delay_s)+1) = v(n);
end

%% implement dynamics filter
x1 = 0; %initialise x1
%read from white noise
for n = 0:N
    x0 = (1-R)*vp(n+1);
    y(n+1) = x0+R*x1;
    x1 = x0;
end

% Slur
tic
yp1 = 0;
for n = N+1:M1
    yp0= y(n-N+1)*C+y(n-N)+yp1*(-C) ;
    y(n+1) = rho*((1-S)*yp0+S*yp1);
    yp1 = yp0;
end

for n = M1+1:M-1
    yp0= y(n-N2+1)*C+y(n-N2)+yp1*(-C) ;
    y(n+1) = rho*((1-S)*yp0+S*yp1);
    yp1 = yp0;
end

%% convolution
IR_guitar = audioread('Taylor814.wav');
y = conv(y,IR_guitar);

%% Normalise
y = y/max(abs(y));
end