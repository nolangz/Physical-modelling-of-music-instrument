%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:Simple harmonic oscillator_Test performance
%%%         Author:Nuolin Lai
%%%         Create Date:04/02/2021
%%%         Last modify date:07/02/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% parameters

SR = 44100;         % sample rate (Hz)
Tf = 1;             % duration (s)
f0 = 500;           % frequency (Hz)
N = 7000;           % Number of frequency 
u0 = ones(1,N);     % initial displacement
v0 = zeros(1,N);    % initial velocity
fv = SR/2*rand(1,N);% frequency vector

alpha = 0.5;  % the scheme free parameter

k = 1/SR;               % time step
w0 = 2*pi*fv;           % angular frequency (rad./s)
Nf = floor(Tf*SR);      % total number of time steps


b = (2-alpha*w0.^2*k^2)./(1 + ((1-alpha)*w0.^2*k^2)/2); %coefficient of u1

% initialize

u2 = u0;                % set initial displacement
u1 = u0+k*v0;           % set second displacement
out = zeros(N,Nf);      % initialse output vector

% main loop

tic
for n=1:Nf
    u = u1.*b - u2;     % scheme update
    out(:,n) = u2;      % write output
    u2 = u1;            % shift state
    u1 = u;             % shift state
end
sum_out = sum(out,1);
toc

%plot
tax = [0:Nf-1]'*k;
plot(tax,sum_out,'k');
xlabel('t');
ylabel('u');
title('Simple Harmonic Oscillator-sum');

%play sound

soundsc(sum_out',SR)
