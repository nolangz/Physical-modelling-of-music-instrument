%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:BTB_Moog_animation_script_sweep_r
%%%         Author:Nuolin Lai
%%%         Create Date:21/02/2021
%%%         Last modify date:021/02/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
clc; 
close all;

%% parameters

h = animatedline;
SR = 44100*4;               % sample rate (Hz)
Tf = 0.05;                   % duration (s)
N =30;
f0 = 2000;       % resonant filter frequency (Hz)
r = linspace(0.1,0.9,N);                    % feedback coeff(choose 0\leq r \leq 1)


%% derived parameter

om0 = 2*pi*f0;         % resonant filter angular frequency(Hz)
Nf  = floor(SR*Tf);    % sample number 
k   = 1/SR;            % time step

%% Initialise

I = eye(4);                                       % identity Matrix
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

z = exp(1i*2*pi*fvec*k);                          % initialise z

He = zeros(Nf,N);                                 % initialise transfer function
Hf_z = zeros(Nf,N);                               % initialise transfer function
Hb_z = zeros(Nf,N);                               % initialise transfer function
Ht_z = zeros(Nf,N);                               % initialise transfer function


%% pre calculate transfer function

for m = 1:N
    A = om0*[-1,0,0,-4*r(m);1,-1,0,0;0,1,-1,0;0,0,1,-1]; % A matrix
    b = om0*[1,0,0,0]';                               % b matrix
    
    %% compute the exact transfer function 
                        
    for n = 1:Nf
        He(n,m) = c'*((2*pi*1i*fvec(n)*I-A)\b);  % Calculate every sample 
    end

    %% compute the Forward transfer function

    for n = 1:Nf
        Hf_z(n,m) = c'/(I*z(n)-(I+k*A))*k*b;      % Calculate every frequency 
    end

    %% compute the Backward transfer function 

    for n = 1:Nf
        Hb_z(n,m) = c'/((I-k*A)*z(n)-I)*k*b*z(n); % Calculate every frequency
    end

    %% compute the Trapezoidal integration transfer function 

    for n = 1:Nf
        Ht_z(n,m) = c'/((I-k*A/2)*z(n)-(I+k*A/2))*(k*b/2*z(n)+k*b/2); % Calculate every frequency
    end

end

%% plot Transfer function in Time domain and Frequency domain
for m = 1:N
    subplot(3,1,1);
    loglog(fvec,abs(Hf_z(:,m)),'r',fvec,abs(He(:,m)),'b','LineWidth',2);
    xlabel('Frequency (Hz)');
    ylabel('magnitude');
    title('Transfer function - Forward vs Exact');
    xlim([20 20000]);
    ylim([10e-6 10e1]);

    subplot(3,1,2);
    loglog(fvec,abs(Hb_z(:,m)),'r',fvec,abs(He(:,m)),'b','LineWidth',2);
    xlabel('Frequency (Hz)');
    ylabel('magnitude');
    title('Transfer function - Backward vs Exact');
    xlim([20 20000]);
    ylim([10e-6 10e1]);

    subplot(3,1,3);
    loglog(fvec,abs(Ht_z(:,m)),'r',fvec,abs(He(:,m)),'b','LineWidth',2);
    xlabel('Frequency (Hz)');
    ylabel('magnitude');
    title('Transfer function - Trapezoid vs Exact');
    xlim([20 20000]);
    ylim([10e-6 10e1]);
    
    drawnow update
end









