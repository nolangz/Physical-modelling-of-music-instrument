%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:MoogVCF_nonlinear_Trapezoidal_Integration
%%%         Author:Nuolin Lai
%%%         Create Date:19/02/2021
%%%         Last modify date:20/02/2021
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

I = eye(4);                                       % identity Matrix
A = om0*[-1,0,0,-4*r;1,-1,0,0;0,1,-1,0;0,0,1,-1]; % A matrix
b = om0*[1,0,0,0]';                               % b matrix
c = [0,0,0,1]';                                   % c matrix

xt = zeros(4,1);                                  % TI state
xtn1 = zeros(4,1);                                % TI next state 
xe = zeros(4,1);                                  % exact state

yt = zeros(Nf,1);                                 % output y from TI

u  = zeros(Nf+1,1)';                              % hold output input sequence
u(2) = 1;                                         % initialise impulse response

tvec = [0:Nf-1]'*k;                               % time vector for plots
fvec = [0:Nf-1]'*SR/Nf;                           % frequency vector for plots

%initialse f
f = [om0*(-tanh(xt(1))+tanh(u(1)-4*r*xt(4)));
    om0*(-tanh(xt(2))+tanh(xt(1)));
    om0*(-tanh(xt(3))+tanh(xt(2)));
    om0*(-tanh(xt(4))+tanh(xt(3)))];


%% main loop to implement Trapezoid method

tic
for n = 1:Nf
    % calculate xtn1 using Newton-Raphason Method
    
    maxIter = 10;   % max number of guesses
    iter = 0;       % iteration counter
    
    % initialise xtn1
    xtn1 = zeros(4,1);
    
    while iter < maxIter
        % calculate fn1
        fn1 = [om0*(-tanh(xtn1(1))+tanh(u(n+1)-4*r*xtn1(4)));
            om0*(-tanh(xtn1(2))+tanh(xtn1(1)));
            om0*(-tanh(xtn1(3))+tanh(xtn1(2)));
            om0*(-tanh(xtn1(4))+tanh(xtn1(3)))];
        % calculate G
        G = xtn1 -xt - k/2*(f+fn1);
        
         % calculate J
        J(1,:) =  [1+k*om0/(2*(cosh(xtn1(1)))^2),0,0,2*k*om0*r/(cosh(u(n)-4*r*xtn1(4)))^2];
        J(2,:) =  [-k*om0/(2*(cosh(xtn1(1)))^2),1+k*om0/(2*(cosh(xtn1(2)))^2),0,0];
        J(3,:) =  [0,-k*om0/(2*(cosh(xtn1(2)))^2),1+k*om0/(2*(cosh(xtn1(3)))^2),0];
        J(4,:) =  [0,0,-k*om0/(2*(cosh(xtn1(3)))^2),1+k*om0/(2*(cosh(xtn1(4)))^2)];
       
        %update xtn1
        xtn1 = xtn1 - J\G;
        iter = iter + 1;
       
    end
    
    yt(n) = c'*xt;                 % write sample to output vector yt 
    
    
    xt = xtn1;                     % shift sample
    f = fn1;
end
simTime = toc;                     % Record simulating time

sound(yt,SR);

%% compute transfer function

Ht = fft(yt);                      % compute transfer function(TI)

%% plot Transfer function in Time domain and Frequency domain
subplot(2,1,1);
plot(tvec,yt,'r','LineWidth',2);
xlabel('Time(s)');
ylabel('magnitude');
title('Transfer function - Time domain');
legend('Trapezoidal Integration')

subplot(2,1,2);
loglog(fvec,abs(Ht),'r','LineWidth',2);
xlabel('Frequency (Hz)');
ylabel('magnitude');
title('Transfer function - Frequency domain');
xlim([20 SR/2]);
ylim([10e-6 10e1]);
legend('Trapezoidal Integration')










