%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:MoogVCF_Transfer Function
%%%         Author:Nuolin Lai
%%%         Create Date:19/02/2021
%%%         Last modify date:19/02/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% parameters

SR = 44100;     % sample rate (Hz)
Tf = 0.2;       % duration (s)
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
xt = zeros(4,1);                                  % TE state
xe = zeros(4,1);                                  % exact state

yf = zeros(Nf,1);                                 % output y from FE
yb = zeros(Nf,1);                                 % output y from BE
yt = zeros(Nf,1);                                 % output y from TE

u  = zeros(Nf,1);                                 % hold output input sequence
u(1) = 1;                                         % initialise impulse response

ut  = zeros(Nf+1,1);                              % hold output input sequence(TE)
ut(1) = 0;                                        % initialise impulse response(TE)
ut(2) = 1;                                        % initialise impulse response(TE)

tvec = [0:Nf-1]'*k;                               % time vector for plots
fvec = [0:Nf-1]'*SR/Nf;                           % frequency vector for plots

z = exp(1i*2*pi*fvec*k);

%% main loop to implement Forward Euler/Forward Euler/Exact Method

tic
for n = 1:Nf
    xfn1  = (I+k*A)*xf+k*b*u(n);   % update state xf from n to n+1 (FE)
    xbn1  = (I-k*A)\(xb+k*b*u(n)); % update state xb from n to n+1 (BE)
    xtn1  = (I-k*A/2)\((I+k*A/2)*xt+k/2*b*(ut(n)+ut(n+1)));   % update state xt from n to n+1 
    yf(n) = c'*xf;                 % write sample to output vector yf (FE)
    yb(n) = c'*xb;                 % write sample to output vector yb (BE)
    yt(n) = c'*xt;                 % write sample to output vector yt (TE)
    xf = xfn1;                     % shift sample(FE)
    xb = xbn1;                     % shift sample(BE)
    xt = xtn1;                     % shift sample(TE)
end
simTime = toc;                     % Record simulating time

%% compute the exact transfer function and 

He = zeros(Nf,1);                          % Initialise Transfer function
for n = 1:Nf
    He(n) = c'*((2*pi*1i*fvec(n)*I-A)\b);  % Calculate every sample 
end
ye = ifft(He,'symmetric');                 % ifft

%% compute transfer function

Hf = fft(yf);                      % compute transfer function(FE)
Hb = fft(yb);                      % compute transfer function(BE)
Ht = fft(yt);                      % compute transfer function(TE)

%% compute the Forward transfer function

Hf_z = zeros(Nf,1);

for n = 1:Nf
    Hf_z(n) = c'/(I*z(n)-(I+k*A))*k*b;      % Calculate every frequency 
end

%% compute the Backward transfer function 

Hb_z = zeros(Nf,1);

for n = 1:Nf
    Hb_z(n) = c'/((I-k*A)*z(n)-I)*k*b*z(n); % Calculate every frequency
end

%% compute the Trapezoidal integration transfer function 

Ht_z = zeros(Nf,1);

for n = 1:Nf
    Ht_z(n) = c'/((I-k*A/2)*z(n)-(I+k*A/2))*(k*b/2*z(n)+k*b/2); % Calculate every frequency
end

%% plot Transfer function in Time domain and Frequency domain
subplot(3,1,1);
loglog(fvec,abs(Hf),'r',fvec,abs(Hf_z),'b','LineWidth',2);
xlabel('Frequency (Hz)');
ylabel('magnitude');
title('Transfer function - Forward');
xlim([20 SR/2]);
ylim([10e-6 10e1]);
legend('numeric method','Z transfer analysis')
hold on;
grid on;

subplot(3,1,2);
loglog(fvec,abs(Hb),'r',fvec,abs(Hb_z),'b','LineWidth',2);
xlabel('Frequency (Hz)');
ylabel('magnitude');
title('Transfer function - Frequency domain');
xlim([20 SR/2]);
ylim([10e-6 10e1]);
legend('numeric method','Z transfer analysis')
hold on;
grid on;

subplot(3,1,3);
loglog(fvec,abs(Ht),'r',fvec,abs(Ht_z),'b','LineWidth',2);
xlabel('Frequency (Hz)');
ylabel('magnitude');
title('Transfer function - Frequency domain');
xlim([20 SR/2]);
ylim([10e-6 10e1]);
legend('numeric method','Z transfer analysis')
hold on;
grid on;

% subplot(3,1,3);
% loglog(fvec,abs(Ht),'r',fvec,abs(H_z),'b','LineWidth',2);
% xlabel('Frequency (Hz)');
% ylabel('magnitude');
% title('Transfer function - Frequency domain');
% xlim([20 SR/2]);
% ylim([10e-6 10e1]);
% legend('numeric method','Z transfer analysis')
% hold on;
% grid on;








