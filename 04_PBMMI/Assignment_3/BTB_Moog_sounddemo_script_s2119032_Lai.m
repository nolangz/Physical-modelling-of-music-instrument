%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:MoogVCF_sound demo_script
%%%         Author:Nuolin Lai
%%%         Create Date:19/02/2021
%%%         Last modify date:20/02/2021
%%%         TRY DIFFERENT INPUT AND PARAMETER, HAVE FUN!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% input parameters
Tf = 5;         % duration (s)
SR = 44100;     % sample rate (Hz)


%% derived parameter

Nf  = floor(SR*Tf);                     % sample number
k   = 1/SR;                             % time step
tvec = [0:Nf-1]'*k;                     % time vector for plots
fvec = [0:Nf-1]'*SR/Nf;                 % frequency vector for plots

%% input signal (choose your input signal) sine/square/triangular/external 
fi = 50;                                 % input frequency
g  = 10;                                 % input gain

%u = sin(2*pi*fi*(0:Nf)/SR);             % sine
%u = sawtooth(2*pi*fi*(0:Nf)/SR);        % triangular 
%u = square(2*pi*fi*(0:Nf)/SR);          % square 
[u,Fs] = audioread('Classic Electric Piano.wav'); % read audio file

Nf = length(u);                        % re-calculate Number of input sample

u = u*g;                                 % implement Gain

%% variable parameter

start_f0 = 100;                         % start f0
end_f0 = 10000;                         % end f0
start_r = 0.9;                          % start r 
end_r =0.4;                             % end r
%f0 = logspace(2,4,Nf+1);                % sweep resonant filter frequency (Hz)
f0 = floor(2100+2000*sawtooth(2*pi*1*(0:Nf)/SR));   
r = linspace(start_r,end_r,Nf+1);       % sweep feedback coeff(choose 0\leq r \leq 1)
om0 = 2*pi.*f0;                         % resonant filter angular frequency(Hz) 

%% Initialise

I = eye(4);                             % identity Matrix
c = [0,0,0,1]';                         % c matrix

xt = zeros(4,1);                        % TI state
xtn1 = zeros(4,1);                      % TI next state 
xe = zeros(4,1);                        % exact state

yt = zeros(Nf,1);                       % output y from TI

tvec = [0:Nf-1]'*k;                     % time vector for plots
fvec = [0:Nf-1]'*SR/Nf;                 % frequency vector for plots

%initialse f
f = [om0(1)*(-tanh(xt(1))+tanh(u(1)-4*r(1)*xt(4)));
    om0(1)*(-tanh(xt(2))+tanh(xt(1)));
    om0(1)*(-tanh(xt(3))+tanh(xt(2)));
    om0(1)*(-tanh(xt(4))+tanh(xt(3)))];

%% main loop to implement Trapezoidal method

tic
for n = 1:Nf
    % calculate xtn1 using Newton-Raphason Method
    
    maxIter = 10;   % max number of guesses
    iter = 0;       % iteration counter
    
    % initialise xtn1
    xtn1 = zeros(4,1);
    
    while iter < maxIter
        % calculate fn1
        fn1 = [om0(n+1)*(-tanh(xtn1(1))+tanh(u(n+1)-4*r(n+1)*xtn1(4)));
            om0(n+1)*(-tanh(xtn1(2))+tanh(xtn1(1)));
            om0(n+1)*(-tanh(xtn1(3))+tanh(xtn1(2)));
            om0(n+1)*(-tanh(xtn1(4))+tanh(xtn1(3)))];
        % calculate G
        G = xtn1 -xt - k/2*(f+fn1);
        % calcualte J
        J = [1+k*om0(n)/(2*(cosh(xtn1(1)))^2),0,0,2*k*om0(n)*r(n)/(cosh(u(n)-4*r(n)*xtn1(4)))^2;
            -k*om0(n)/(2*(cosh(xtn1(1)))^2),1+k*om0(n)/(2*(cosh(xtn1(2)))^2),0,0;
            0,-k*om0(n)/(2*(cosh(xtn1(2)))^2),1+k*om0(n)/(2*(cosh(xtn1(3)))^2),0
            0,0,-k*om0(n)/(2*(cosh(xtn1(3)))^2),1+k*om0(n)/(2*(cosh(xtn1(4)))^2)
            ];
       
        %update
       xtn1 = xtn1 - J\G;
       iter = iter + 1;
       
    end
    
    yt(n) = c'*xt;                 % write sample to output vector yt 
    
    xt = xtn1;                     % shift sample
    f = fn1;
end
simTime = toc;                     % Record simulating time


%% compute transfer function

Ht = fft(yt);                      % compute transfer function(TI)

%% plot Transfer function in Time domain and Frequency domain
subplot(2,1,1);
plot(tvec,yt,'r','LineWidth',1);
xlabel('Time(s)');
ylabel('magnitude');
title('Transfer function - Time domain');
legend('Non-Linear Trapezoidal Integration')

subplot(2,1,2);
loglog(fvec,abs(Ht),'r','LineWidth',1);
xlabel('Frequency (Hz)');
ylabel('magnitude');
title('Transfer function - Frequency domain');
xlim([20 SR/2]);
legend('Non-Linear Trapezoidal Integration')

%% Normalise and play
yt = yt/max(abs(yt));
sound(yt,SR);










