%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:Simple harmonic oscillator_conserved energy
%%%         Author:Nuolin Lai
%%%         Create Date:08/02/2021
%%%         Last modify date:08/02/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% parameters

SR = 44100;     % sample rate (Hz)
Tf = 1;         % duration (s)
f0 = 15e3;      % frequency (Hz)
u0 = 1;         % initial displacement
v0 = 0;         % initial velocity


k = 1/SR;               % time step
w0 = 2*pi*f0;           % angular frequency (rad./s)
Nf = floor(Tf*SR);      % total number of time steps

alpha = 2/w0^2/k^2-cos(w0*k)/(1-cos(w0*k));  % the scheme free parameter

b = (2-alpha*w0^2*k^2)/(1 + ((1-alpha)*w0^2*k^2)/2); %coefficient of u1

be = 2*cos(w0*k);                                    % coefficient of u1（exact）

% initialize

%assert
if(k > 2/(w0*sqrt(2*alpha-1)) && alpha >= 0.5)
    error('Stability condition violated')
end

u2 = u0;                        % set initial displacement
u1 = u0+k*v0+k^2/2*(-w0^2*u0);  % set second displacement

u2e = u0;                       % set initial displacement
u1e = u0+k*v0+k^2/2*(-w0^2*u0); % set second displacement


out = zeros(Nf,1);      % output vector
outex = zeros(Nf,1);    % output vector（exact）
er = zeros(Nf,1);       % error vector
he = zeros(Nf,1);        % energy vector
% main loop

tic
for n=1:Nf
    u = b*u1 - u2;        % scheme update
    out(n) = u2;          % write output
    u2 = u1;              % shift state
    u1 = u;
    ue = be*u1e - u2e;    % scheme update
    outex(n) = u2e;       % write output
    he(n) = 1/(2*k^2)*(u1e^2+(w0^2*k^2-2)*u1e*u2e+u2e^2);
    u2e = u1e;            % shift state
    u1e = ue;
    er(n) = out(n)-outex(n);
end
toc

% plot

tax = [0:Nf-1]'*k;
plot(tax,he,'k');
xlabel('t');
ylabel('u');
title('Simple Harmonic Oscillator-energy');

% play sound

soundsc(out,SR)
