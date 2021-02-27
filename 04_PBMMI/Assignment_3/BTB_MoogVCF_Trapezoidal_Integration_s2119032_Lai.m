%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:BTB_MoogVCF_Trapezoidal_Integration
%%%         Author:Nuolin Lai
%%%         Create Date:04/02/2021
%%%         Last modify date:04/02/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% parameters

SR = 44100;     % sample rate (Hz)
Tf = 0.05;       % duration (s)
f0 = 1000;       % resonant filter frequency (Hz)
r = 0.7;        % feedback coeff(choose 0\leq r \leq 1)

%% derived parameter

om0 = 2*pi*f0;         % resonant filter angular frequency(Hz)
Nf  = floor(SR*Tf);    % sample number 
k   = 1/SR;            % time step
a   = sqrt(2)*r^(1/4); % magnitude coefficient

%% Initialise
A = om0*[-1,0,0,-4*r;1,-1,0,0;0,1,-1,0;0,0,1,-1]; % A matrix
I = eye(4);                                       % identity Matrix
b = om0*[1,0,0,0]';                               % b matrix
c = [0,0,0,1]';                                   % c matrix

xt = zeros(4,1);                                  % TI state
yt = zeros(Nf,1);                                 % output y from TI

u  = zeros(Nf+1,1);                               % hold output input sequence
u(2) = 1;                                         % initialise impulse response

tvec = [0:Nf-1]'*k;                               % time vector for plots
fvec = [0:Nf-1]'*SR/Nf;                           % frequency vector for plots

%% main loop to implement Trapezoidal Method

tic
for n = 1:Nf
    xtn1  = (I-k*A/2)\((I+k*A/2)*xt+k/2*b*(u(n)+u(n+1)));   % update state xt from n to n+1 
    yt(n) = c'*xt;                 % write sample to output vector yt 
    xt = xtn1;                     % shift sample(BE)
end
simTime = toc;                     % Record simulating time

%% compute the exact transfer function and ifft

He = zeros(Nf,1);                          % Initialise Transfer function
for n = 1:Nf
    He(n) = c'*((2*pi*1i*fvec(n)*I-A)\b);  % Calculate every sample 
end
ye = ifft(He,'symmetric');                 % ifft

%% compute transfer function

Ht = fft(yt);                      % compute transfer function(TI)

%% plot Transfer function in Time domain and Frequency domain
subplot(2,1,1);
plot(tvec,yt,'r',tvec,ye,'k','LineWidth',2);
xlabel('Time(s)');
ylabel('magnitude');
title('Transfer function - Time domain');
legend('Trapezoidal Integration','Exact')

subplot(2,1,2);
loglog(fvec,abs(Ht),'r',fvec,abs(He),'k','LineWidth',2);
xlabel('Frequency (Hz)');
ylabel('magnitude');
title('Transfer function - Frequency domain');
xlim([20 SR/2]);
ylim([10e-6 10e1]);
legend('Trapezoidal Integration','Exact')










