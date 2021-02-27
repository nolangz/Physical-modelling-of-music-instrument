%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:Simple harmonic oscillator
%%%         Author:Nuolin Lai
%%%         Create Date:04/02/2021
%%%         Last modify date:04/02/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% parameters

SR = 44100;     % sample rate (Hz)
Tf = 0.05;      % duration (s)
f0 = 1000;      % resonant filter frequency (Hz)
r = 0.7;        % feedback coeff(choose 0\leq r \leq 1)
alpha = 0.9;    % weighted parameter

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

xp = zeros(4,1);                                  % PI state
yp = zeros(Nf,1);                                 % output y from PI

u  = zeros(Nf+1,1);                               % hold output input sequence
u(2) = 1;                                         % initialise impulse response

tvec = [0:Nf-1]'*k;                               % time vector for plots
fvec = [0:Nf-1]'*SR/Nf;                           % frequency vector for plots

%% main loop to implement Parameterised method

tic
for n = 1:Nf
    xpn1  = (I-alpha*k*A/2)\((I+k*A-k/2*alpha*A)*xp+k*b/2*(u(n)+u(n+1)));   % update state xp from n to n+1 
    yp(n) = c'*xp;                 % write sample to output vector ytp
    xp = xpn1;                     % shift sample(PE)
end
simTime = toc;                     % Record simulating time

%% compute the exact transfer function 

He = zeros(Nf,1);                          % Initialise Transfer function
for n = 1:Nf
    He(n) = c'*((2*pi*1i*fvec(n)*I-A)\b);  % Calculate every sample 
end
ye = ifft(He,'symmetric');

%% compute transfer function

Hp = fft(yp);                      % compute transfer function(PI)

%% plot Transfer function in Time domain and Frequency domain
subplot(2,1,1);
plot(tvec,yp,'r',tvec,ye,'k','LineWidth',2);
xlabel('Time(s)');
ylabel('magnitude');
title('Transfer function - Time domain');
legend('Parameterised','Exact')

subplot(2,1,2);
loglog(fvec,abs(Hp),'r',fvec,abs(He),'k','LineWidth',2);
xlabel('Frequency (Hz)');
ylabel('magnitude');
title('Transfer function - Frequency domain');
xlim([20 SR/2]);
ylim([10e-6 10e1]);
legend('Parameterised','Exact')










