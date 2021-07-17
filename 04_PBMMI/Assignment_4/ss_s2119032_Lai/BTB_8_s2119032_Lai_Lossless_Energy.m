%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:Stiff string Frequency dependent loss
%%%         Author:Nuolin Lai
%%%         Create Date:11/03/2021
%%%         Last modify date:11/03/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%%%%% flags

plot_on = 0;                % in-loop plotting on (1) or off (0)
itype = 1;                  % type of input: 1: struck, 2: plucked

%%%%% parameters

% physical string parameters

T = 60;                     % tension (N)
r = 0.0004;                 % string radius (m)
E = 2e11;                   % Young's modulus (Pa) for steel
rho = 7850;                 % density (kg/m^3)
T60 = 5;                    % T60 (s)
L = 1;                      % length (m)

% I/O

SR = 44100;                 % sample rate (Hz)
Tf = 5;                     % duration of simulation (s)

xi = 0.9;                   % coordinate of excitation (normalised, 0-1)
famp = 1;                   % peak amplitude of excitation (N)
dur = 0.001;                % duration of excitation (s)
exc_st = 0.1;               % start time of excitation (s)

xo = 0.1;                   % coordinate of output (normalised, 0-1)

%%%%% Q1: perform standard checking on all these parameters (for non-negativity, e.g.), producing an
%%%%% exit error message or warning as appropriate. In addition, find a way
%%%%% to check that, from these parameters, the number N of segments which
%%%%% will be later derived is less than 10000, and exit if it is violated. 
%%%%% Also, find a way to check that, for a given choice of xo or xi, the
%%%%% resulting grid location, in metres, will be at least h metres away 
%%%%% from either endpoint of the string, at x=0 or x=L. Also check that
%%%%% the full duration of your excitation pulse falls within the time span
%%%%% [0, Tf] of the simulation!

%%%%% Add this code below:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% derived parameters

A = pi*r^2;                 % string cross-sectional area
I = 0.25*pi*r^4;            % string moment of inertia
c = sqrt(T/(rho*A));        % wave speed
K = sqrt(E*I/(rho*A));      % stiffness constant 
sig = 6*log(10)/T60;        % frequency independent loss parameter

%%%%% grid

k = 1/SR;                   % time step

%%%%% Q2: determine, from the system parameters c and K and the time step
%%%%% k, a minimal grid spacing hmin. Then, from hmin and L, determine the
%%%%% number of grid spacings N into which the string may be divided evenly,
%%%%% as close to hmin (but not below!) as possible. Then, reset h from N and
%%%%% L. 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% YOUR CODE (SEE Q2)

hmin = sqrt(1/2*(c^2*k^2+sqrt(c^4*k^4+16*K^2*k^2)));
N = floor(L/(hmin));
h = L/N;

%% stability assert

assert(N<10000,'error! Please check your input parameter');
assert(h >= hmin,'error! Please check your input parameter');
assert(sig>=0,'error! Please check your input parameter');
assert((xi*L)>h),'error! Please check your input parameter';
assert((L-xi*L)>h,'error! Please check your input parameter');
assert(dur<Tf,'error! Please check your input parameter');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create force signal

Nf = floor(SR*Tf);                      % number of time steps

%%%%% Q3: from the itype flag, and the parameters dur, exc_st and famp,
%%%%% generate a vector f of length Nf samples, representing the input
%%%%% force in Newtons. 

%%%%% Add this code below:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_vec = (0:Nf-1)/SR;
f = zeros(1,Nf);

for n = 1:Nf
    if t_vec(n)>=exc_st && t_vec(n)<=exc_st+dur
        f(n) = 1/2*famp*(1-cos(itype*pi()*(t_vec(n)-exc_st)/dur));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% scheme coefficients

% matrices

%%%%% Q4: Referring to the matrix form of the update in the main loop, of
%%%%% the form u = B*u1-C*u2+J*f(n), design (N-1)x(N-1) matrices B and C,
%%%%% in sparse form, as well as an (N-1)x1 vector J, selecting the input
%%%%% location. Also, for the output y(n) = c'*u, define an (N-1)x1 vector
%%%%% c, selecting the output location. All of these matrices and vectors should be
%%%%% represented in sparse form!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% YOUR CODE (SEE Q4)
I_mtx = speye(N-1);
e = ones(N-1,1);

Dx = 1/h * spdiags([e -e],0:1,N-1,N-1); % calculate Dx for main loop

Dxx = 1/h^2*spdiags([e -2*e e], -1:1, N-1,N-1);
spdiags_1 = spdiags([e -4*e 6*e -4*e e], -2:2, N-1,N-1);
spdiags_1(1,1) = 5;
spdiags_1(N-1,N-1) = 5;
Dxxxx = (1/h^4)*spdiags_1;
B = 2 * I_mtx+c^2*k^2*Dxx;
C = I_mtx;
li = floor(xi/h);
lo = floor(xo/h);
J = zeros(N-1,1);
J(li) = 1;
J = k^2/(rho*A*h)*J;
c_v = zeros(N-1,1);
c_v(lo) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% initialise scheme variables

u2 = zeros(N-1,1);              % state
u1 = u2;                        % state
u = u2;                         % state

y = zeros(Nf,1);                % output
H = zeros(Nf,1);                % energy

%%%%% main loop

tic
for n=1:Nf
    
    % update state, and insert current value of f.
    
    u = B*u1-C*u2+J*f(n);
    
    H(n) = 1/2*norm((u -u1)/k)^2 + c^2/2*norm(Dx*u)^2 +K^2/2*norm(Dxx*u)^2;
    % read output
    
    y(n) = c_v'*u;
    
    % plot
    
    if(plot_on==1)
        % draw current state
        figure(1)
        plot([1:N-1]'*h, u, 'k');
        xlabel('x (m)')
        ylabel('u (m)')
        axis([0 L -0.005 0.005])
        drawnow
    end
    
    % shift state
    
    u2 = u1;
    u1 = u;
    
end
toc

%%%%% play sound

soundsc(y,SR);

%%%%% plot spectrum

figure(2)
subplot(2,1,1)
yfft = 10*log10(abs(fft(y)));
plot([0:Nf-1]'*(SR/Nf), yfft, 'k')
xlim([0 SR/2])
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Transform of output')
subplot(2,1,2)
plot(t_vec,H)

%%%%% BEYOND THE BASICS
%%%%% There are many ways to go from this point. You might consider:

%%%%% B1: find physical string parameters (see above) corresponding to an acoustic guitar, with EADGBE tuning. 
%%%%% You may do this partly by tracking down published values in the literature. 
%%%%% Length L, Young's modulus E, radius r and density rho will be easy to find. Tension T will be more difficult. 
%%%%% T60 will be nearly impossible to find. For any such values which you cannot find, find recordings of
%%%%% isolated notes, and attempt to estimate these (scientifically, using Fourier transform analysis to ascertain the location of spectral
%%%%% peaks and bandwidths). You will need to find six separate sets of
%%%%% parameters for each string on the guitar (EADGBE). Add code to toggle
%%%%% the different notes (say, a flag string_num taking on values from 1
%%%%% to 6) in order to select a particular parameter set. 

%%%%% B2: Generalise this code such that it is capable of emulating not one string, 
%%%%% but all six in the acoustic guitar. It should produce a monophonic
%%%%% audio output, corresponding to a direct sum of the individual output signals from
%%%%% each string. Please configure your code such that it produces a
%%%%% single vector (which you can still call "y"), of duration 5 seconds,
%%%%% which consists of the six notes EADGBG played in ascending order, and
%%%%% at intervals of 0.5 seconds, with the first note occuring at t=0 seconds. You may decide whether a) you simply
%%%%% run six main loops serially (i.e., one after another), generating six output vectors 
%%%%% which you then simply add together or b) whether
%%%%% you will use one main loop, over which you perform the updates on all
%%%%% strings. NB: you will need to add separate excitation parameters,
%%%%% i.e., xi, famp, exc_st and dur, as well as output location xo for
%%%%% each separate string! Can you consolidate the state of all six
%%%%% strings into a single vector? This allows a single state space update
%%%%% of the entire system.

%%%%% B3: if you are working on a six (or more) string...consider ways of mixing the six output vectors from the six
%%%%% strings to a stereo output (you could imagine panning these channels
%%%%% individually from -1 (left) to 1 (right))

%%%%% B4: if you are working on a six (or more) string, create a way of managing or generating many input events (beyond just simply
%%%%% playing six notes in succession). Do not forget that each note will
%%%%% be described by five parameters: xi, famp, dur and exc_st, as well as
%%%%% the string number (1-6 in this case). Can you develop functions which
%%%%% emulate upward or downward strums? A block chord? Can you employ
%%%%% targeted randomness in your note parameter generation to get more natural
%%%%% sound?

%%%%% B5: find a more perceptually reasonable way of setting string
%%%%% parameters (e.g., in terms of fundamental frequency and
%%%%% inharmonicity, or MIDI note number, or even note names like A or B or material names like steel or gold, referring to tables where material values are held)

%%%%% B6: generate the most interesting piece of music (of exactly 20
%%%%% seconds duration) that you can. Best examples in descending order will receive 3,2,and 1 extra points,
%%%%% respectively. Include this as a separate .wav file in your
%%%%% submission...if you dare! If it's really bad, I reserve the right to
%%%%% give negative points. For example...it had better at least be in
%%%%% stereo. 

%%%%% B7: introduce frequency-dependent loss (as well as a way of
%%%%% specifying a variation in T60). This is covered in NSS. If you do
%%%%% this, use an explicit scheme, and make sure that you adjust the
%%%%% numerical stability condition, from Q2, accordingly! Your C and B matrices, from Q4, will
%%%%% also need minor changes!

%%%%% B8: Build in a calculation of the numerical energy for the scheme,
%%%%% following the definitions in lecture and NSS---NB the energy will be
%%%%% conserved to machine accuracy only under lossless conditions! You
%%%%% could plot the relative variation in the total energy as a function
%%%%% of time step, and demonstrate that it is conserved to machine
%%%%% precision. 





