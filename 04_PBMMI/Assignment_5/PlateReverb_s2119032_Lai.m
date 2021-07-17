%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:PLATE REVERB
%%%         Author:Nuolin Lai
%%%         Create Date:21/03/2021
%%%         Last modify date:24/03/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% Gobal parameters

SR = 44100;     % sample rate (Hz)
Tf = 5;         % duration (s)
r = 0.7;        % feedback coeff(choose 0\leq r \leq 1)
Nf = SR * Tf;   % number of sample for impulse response

flag = 1;       %（0：default;1:Aluminium）


% F = zeros(Nf,1);
% F(1) = 1;

[F, Fs] = audioread('PlateReverb_s2119032_Lai_dry.wav');


%% physical parameter


switch flag
    case 0
    Lx = 2;               % Length of the plate
    Ly = 1;               % Width of the plate
    T = 700;              % Tension
    H = 0.5e-3;           % thickness
    rho = 7870;           % density
    E = 2e11;             % young's modular
    v = 0.3;              % Poisson's ratio
    T60_min = 1;          % minimal T60
    T60_max = 4;          % maximal T60
    case 1  %Aluminium
    Lx = 3;               % Length of the plate
    Ly = 2;               % Width of the plate
    T = 700;              % Tension
    H = 0.8e-3;           % thickness
    rho = 2700;           % density
    E = 7e11;             % young's modular
    v = 0.33;             % Poisson's ratio
    T60_min = 1.5;        % minimal T60
    T60_max = 6;          % maximal T60
end

        

%% derived parameter

K = sqrt(E*H^2/(12*rho*(1-v^2)));   % resonant filter angular frequency(Hz)
c = sqrt(T/(rho*H));                % wave speed
sigma_min = 6*log(10)/T60_max;      % loss parameter minimal
sigma_max = 6*log(10)/T60_min;      % loss parameter maximal
k   = 1/SR;                         % time step

%% stability condition

omega_max = 2/k;                    % stability condition
beta_max = sqrt((-c^2 + sqrt(c^4 + 4*K^2 * omega_max^2))/(2*K^2)); % max wave number
q_min = [1,1];                       
beta_min = sqrt((q_min(1)*pi()/Lx)^2 + (q_min(2)*pi()/Ly)^2);      % min wave number
e1 = (sigma_max - sigma_min)/((beta_max)^2-(beta_min)^2);          % loss constant
e0 =  sigma_min - e1 *beta_min^2;                                  % loss constant

%% meshgrid
Qx = floor(sqrt(beta_max^2-(pi/Ly)^2)* Lx/pi());    % number of grid in X
Qy = floor(sqrt(beta_max^2-(pi/Lx)^2)* Ly/pi());    % number of grid in Y

Q = Qx *Qy;                                         % total number of meshgrid

[x,y] = meshgrid(1:Qx,1:Qy);                        % create meshgrid

qx = reshape(x,[Q,1]);                              % two colume vectors
qy = reshape(y,[Q,1]);                              % two colume vectors

beta_q = sqrt((qx*pi()/Lx).^2+(qy*pi()/Ly).^2);     % wave number according to q
pass = beta_q >= 0 & beta_q < beta_max;
beta_q = beta_q(pass);                              % select the pass term for stability

%% initial I/O

input = [1,0.5];     % in
outputL = [0.9,0.5]; % left out
outputR = [1.5,0.2]; % right out

% modal shape
Phi = 2*sin(qx*pi()*input(1)/Lx).*sin(qy*pi()*input(2)/Ly)/sqrt(Lx*Ly);
Phi_L = 2*sin(qx*pi()*outputL(1)/Lx).*sin(qy*pi()*outputL(2)/Ly)/sqrt(Lx*Ly);
Phi_R = 2*sin(qx*pi()*outputR(1)/Lx).*sin(qy*pi()*outputR(2)/Ly)/sqrt(Lx*Ly);

% select the pass term for stability
Phi = Phi(pass);
Phi_L = Phi_L(pass);
Phi_R = Phi_R(pass);

% new mode number
Np = length(Phi);

%% Finite difference scheme parameter calculation
OMEGA = sqrt(c^2 * beta_q.^2 + K^2 * beta_q.^4);
SIGMA = e0 + e1 * beta_q.^2;


% coefficient vector
a = (2 - k^2 * OMEGA.^2) ./( 1 + k * SIGMA);
b = (k * SIGMA -1 ) ./(k * SIGMA +1);
c = k^2 ./(1 + SIGMA * k);

p1 = zeros(Np,1);
p2 = zeros(Np,1);

N = length(F);

u = zeros(N,2);

%% main loop to implement Finite diference scheme

tic
for n = 1:N
    
    p = a .* p1 + b.* p2 + c .* Phi * F(n); 
    u(n,1) = Phi_L' * p;
    u(n,2) = Phi_R' * p;
    p2 = p1;
    p1 = p;
    
end
simTime = toc;                     % Record simulating time

%% Normalisation
max_u = max(max(abs(u)));
u = u/max_u;

%% Play

sound(u,SR);
audiowrite('output.wav',u,SR);

%% spectrogram
u = sum(u,2);
[mag, fre , t] = spectrogram(u,4092,128,4092,SR);
t_vec = (0:N-1)/SR;

subplot(2,1,1);
plot(t_vec,u);
title('Waveform')
xlabel('time(s)');
ylabel('magnitude');

subplot(2,1,2);
h = pcolor(t,fre,20*log10((abs(mag))));
h.EdgeColor = 'none';
title('Spectrogram')
xlabel('time(s)');
ylabel('frequency(Hz)');
colorbar;
set(gca,'Layer','top');
set(gca,'YScale','log','YDir','normal','YTick',[20,50,100,200,500,1000,2000,5000,10000,20000]);
axis tight 

figure(2)
h = pcolor(t,fre,20*log10((abs(mag))));
h.EdgeColor = 'none';
title('Spectrogram')
xlabel('time(s)');
ylabel('frequency(Hz)');
colorbar;
set(gca,'Layer','top');
set(gca,'YScale','log','YDir','normal','YTick',[20,50,100,200,500,1000,2000,5000,10000,20000]);
axis tight 