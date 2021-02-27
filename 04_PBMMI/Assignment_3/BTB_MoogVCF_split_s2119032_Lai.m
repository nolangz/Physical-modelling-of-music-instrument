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
Tf = 0.05;       % duration (s)
f0 = 1000;       % resonant filter frequency (Hz)
r = 0.5;        % feedback coeff(choose 0\leq r \leq 1)

%% derived parameter

om0 = 2*pi*f0;         % resonant filter angular frequency(Hz)
Nf  = floor(SR*Tf);    % sample number 
k   = 1/SR;            % time step
a   = sqrt(2)*r^(1/4); % magnitude coefficient

%% Initialise
A = om0*[-1,0,0,-4*r;1,-1,0,0;0,1,-1,0;0,0,1,-1];  % A matrix
A2 = om0*[-rand(1),0,0,-4*rand(1);rand(1),-rand(1),0,0;0,rand(1),-rand(1),0;0,0,rand(1),-rand(1)];% A2 matrix
A1 = A-A2;                                         % A1 matrix
I = eye(4);                                        % identity Matrix
b = om0*[1,0,0,0]';                                % b matrix
c = [0,0,0,1]';                                    % c matrix

xs = zeros(4,1);                                   % split state
ys = zeros(Nf,1);                                  % output y from split

u  = zeros(Nf+1,1);                                % hold output input sequence
u(2) = 1;                                          % initialise impulse response

tvec = [0:Nf-1]'*k;                                % time vector for plots
fvec = [0:Nf-1]'*SR/Nf;                            % frequency vector for plots

%% main loop to implement split Method

tic
for n = 1:Nf
    xsn1  = (I-k*A1/2)\((I+k*A1/2+k*A2)*xs+k/2*b*(u(n)+u(n+1)));   % update state xs from n to n+1 
    ys(n) = c'*xs;                 % write sample to output vector ys 
    xs = xsn1;                     % shift sample(BE)
end
simTime = toc;                     % Record simulating time

%% compute the exact transfer function and 

He = zeros(Nf,1);                          % Initialise Transfer function
for n = 1:Nf
    He(n) = c'*((2*pi*1i*fvec(n)*I-A)\b);  % Calculate every sample 
end
ye = ifft(He,'symmetric');                 % ifft

%% compute transfer function

Hs = fft(ys);                      % compute transfer function(SI)

%% plot Transfer function in Time domain and Frequency domain
subplot(2,1,1);
plot(tvec,ys,'r',tvec,ye,'k','LineWidth',2);
xlabel('Time(s)');
ylabel('magnitude');
title('Transfer function - Time domain');
legend('spilt','Exact')

subplot(2,1,2);
loglog(fvec,abs(Hs),'r',fvec,abs(He),'k','LineWidth',2);
xlabel('Frequency (Hz)');
ylabel('magnitude');
title('Transfer function - Frequency domain');
xlim([20 SR/2]);
ylim([10e-6 10e1]);
legend('spilt','Exact')










