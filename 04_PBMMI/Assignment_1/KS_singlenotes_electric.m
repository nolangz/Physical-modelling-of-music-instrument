%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:KS_singlenotes_electric
%%%         Author:Nuolin Lai
%%%         Create Date:28/01/2020
%%%         Last modify date:18/01/2020
%%%         Input BPM,syncopated note type(1,1/2,1/4),fundamental freqeuncy
%%%         Sampling frequency,Dynamics parameter,decay stretching factor,
%%%         pick position
%%%         BPM = 60~120 (suggestion)
%%%         syn1 = 4,2,1,1/2,1/4,1/8....
%%%         f0   = frequency of notes
%%%         Fs: sampling rate
%%%         R（0，1）: dynamics parameter R
%%%         S（0，1）: the decay stretching factor
%%%         Mu (0,1): pick position (0,1) 1/2 is in the middle
%%%         example KS_singlenotes_electric(70,3,110,44100,0.95,0.9,0.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function y = KS_singlenotes_eletric(BPM,syn,f0,Fs,R,S,Mu);
%set global paramenter
Bt  = 60/BPM;          % beat duration in seconds.
factor = 1;            % correct factor between feeling duration and T60
T60 = Bt*syn*factor;   % duration for syncopated note
Fs  = 44.1e3;          % sampling rate in Hz
tau = T60/log(1000);   % calculate tau         
rho =  exp(-1/(f0*tau))/abs(cos(pi*f0/Fs)); %loss parameter rho
P1  = Fs/f0;           % the real value for the period of the first partial
Pa = S;
%omega0 = 2*pi*f0; 
% if f0 <= Fs/40
%     Pa = S;
% else
%     Pa = 1/(tan(-S*sin(omega0/Fs)/((1-S)-S*cos(omega0/Fs)))*omega0/Fs);
% end
Nexact = P1-Pa;             % ideal number of samples
M = floor(T60*Fs);          % the length of output 
N = floor(Nexact);          % number of sample in delay line
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
%% feedback and distortion parameter

delay_fb   = 0.01;                  %feedback delay from speaker in second
delay_fbs  = delay_fb*Fs;    %feedback delay from speaker in sample number
pre_gain   = 10;                     %pre distortion gain 
dist_gain  = 0.8;                   %distortion output level
clean_gain = 1-dist_gain;           %pre-distortion output level
fb_gain    = 0.0005;                   %feedback gain 
fb         = zeros(M+delay_fbs,1);  %initialise feedback vector
y_out      = zeros(M,1);            %initialise output vector
%% distortion function
    function y = softclipping(x);
            if x >= 1
                y = 2/3;
            elseif x <= -1
                y = -2/3;
            else
                y = x+x^3/3;
            end
    end
%% implement dynamics filter

x1 = 0; %initialise x1

%read from white noise
for n = 0:N
    x0 = (1-R)*vp(n+1);
    y(n+1) = x0+R*x1;
    x1 = x0;
end

%decay stretching and tuning with feedback and distortion
yp1 = 0;
for n = N+1:M-1
    yp0= y(n-N+1)*C+y(n-N)+yp1*(-C);
    y(n+1) = rho*((1-S)*yp0+S*yp1)+fb(n+1);
    yp1 = yp0;
    %feedback and distortion
    y_dist = softclipping(pre_gain*y(n+1));
    y_out(n+1)  = y(n+1)*clean_gain + y_dist*dist_gain;
    fb(n+delay_fbs+1)     = y_dist*fb_gain;
end

%% Normalise
y = y_out/max(abs(y_out));
end