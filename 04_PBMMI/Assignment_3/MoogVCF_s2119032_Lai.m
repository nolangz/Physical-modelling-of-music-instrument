%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:MoogVCF_Basic
%%%         Author:Nuolin Lai
%%%         Create Date:18/02/2021
%%%         Last modify date:18/02/2021
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

%% assert

%take the square of scheme 26 and solve for k

for m = 1:4
    assert(k <= 2*(1-a*cos(pi/4+pi/2*m))/(om0*(a^2-2*a*cos(pi/4+pi/2*m)+1)));
end

%% Initialise

I = eye(4);                                       % identity Matrix
A = om0*[-1,0,0,-4*r;1,-1,0,0;0,1,-1,0;0,0,1,-1]; % A matrix
b = om0*[1,0,0,0]';                               % b matrix
c = [0,0,0,1]';                                   % c matrix

xf = zeros(4,1);                                  % FE state
xb = zeros(4,1);                                  % BE state
xe = zeros(4,1);                                  % exact state

yf = zeros(Nf,1);                                 % output y from FE
yb = zeros(Nf,1);                                 % output y from BE

u  = zeros(Nf+1,1);                                 % hold output input sequence
u(1) = 1;                                         % initialise impulse response

tvec = [0:Nf-1]'*k;                               % time vector for plots
fvec = [0:Nf-1]'*SR/Nf;                           % frequency vector for plots

%% main loop to implement Forward Euler/Forward Euler/Exact Method

tic
for n = 1:Nf
    xfn1  = (I+k*A)*xf+k*b*u(n);   % update state xf from n to n+1 (FE)
    xbn1  = (I-k*A)\(xb+k*b*u(n)); % update state xb from n to n+1 (BE)
    yf(n) = c'*xf;                 % write sample to output vector yf (FE)
    yb(n) = c'*xb;                 % write sample to output vector yb (BE)
    xf = xfn1;                     % shift sample(FE)
    xb = xbn1;                     % shift sample(BE)
end
simTime = toc;                     % Record simulating time

%% compute the exact transfer function and 

He = zeros(Nf,1);                          % Initialise Transfer function
for n = 1:Nf
    He(n) = c'*((2*pi*1i*fvec(n)*I-A)\b);  % Calculate every sample 
end
ye = ifft(He,'symmetric');

%% compute transfer function

Hf = fft(yf);                      % compute transfer function(FE)
Hb = fft(yb);                      % compute transfer function(BE)

%% plot Transfer function in Time domain and Frequency domain
subplot(2,1,1);
p = plot(tvec,yf,'r',tvec,yb,'b',tvec,ye,'k','LineWidth',2);
xlabel('Time(s)');
ylabel('magnitude');
title('Transfer function - Time domain');
legend('forward Euler','backward Eular','Exact')

subplot(2,1,2);
loglog(fvec,abs(Hf),'r',fvec,abs(Hb),'b',fvec,abs(He),'k','LineWidth',2);
xlabel('Frequency (Hz)');
ylabel('magnitude');
title('Transfer function - Frequency domain');
xlim([20 SR/2]);
ylim([10e-6 10e1]);
legend('forward Euler','backward Eular','Exact')










