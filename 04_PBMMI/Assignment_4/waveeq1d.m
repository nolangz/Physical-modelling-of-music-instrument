% 1D wave equation
% S. Bilbao, 1 December 2020

clear all
close all

% flags

plot_on = 0;        % plot if = 1

% parameters

SR = 44100;         % sample rate (Hz)
L = 1;              % string length (m)
T = 40;             % tension (N)
rho = 7850;         % density (kg/m^3)...for steel
r = 0.0005;         % string radius (m)...typical of musical strings
Tf = 5;             % simulation duration (s)
xi = 0.7;           % center of initial condition (m)
wid = 0.05;         % width of initial condition (m)
u0amp = 0.01;       % amplitude of initial condition (m)
lambda = 1;         % nominal Courant number

% derived parameters

A = pi*r^2;         % string cross sectional area (m^2)
c = sqrt(T/(rho*A));% wave speed (m/s)
k = 1/SR;           % time step (s)

h = c*k/lambda;     % grid spacing
N = floor(L/h);     % set integer number of segments
h = L/N;            % reset grid spacing  
lambda = c*k/h;     % reset Courant number

Nf = floor(Tf*SR);  % number of time steps

% create second derivative operator
e = ones(N-1,1);
Dxx = (1/h^2)*spdiags([e -2*e e], -1:1, N-1,N-1);
B = 2*speye(N-1)+c^2*k^2*Dxx;


% initial condition setup

u0 = zeros(N-1,1);
for qq=1:N-1
    x = (qq)*h;
    dist = abs(x-xi);
    if(dist<=wid)
        u0(qq) = 0.5*u0amp*(1+cos(pi*dist/wid));
    end
end

% initialize

u1 = u0;
u2 = u0;
u = zeros(N-1,1);

% main loop
x = zeros(Nf, 1);
for n=1:Nf
    
    u = B*u1-u2;
    
    % plot 
    
    if(plot_on==1)
        plot([1:N-1]*h, u, 'k')
        axis([0 L -u0amp u0amp])
        xlabel('x')
        ylabel('u')
        title('Ideal String')
        drawnow
    end
    u2 = u1;
    u1 = u;
    x(n) = u(100);
end


