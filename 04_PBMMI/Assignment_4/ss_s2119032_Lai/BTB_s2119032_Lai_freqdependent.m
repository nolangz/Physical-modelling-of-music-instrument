%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:Stiff string Frequency dependent loss function
%%%         Author:Nuolin Lai
%%%         Create Date:13/03/2021
%%%         Last modify date:13/03/2021
%%%         Include:
%%%         1. EADGBE tuning
%%%         2. switch string material and radius
%%%         3. input string number/Fret number/stimulating
%%%         parameter/duration/inharmonity parameter to get a corresponding
%%%         sound
%%%         4.add randomness to improve natureness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y  =  BTB_s2119032_Lai_freqdependent(xi,xo,famp,dur,exc_st,string_n,itype,r,material,B,position,Tf)

SR = 44100;                                  % sample rate (Hz)
Nf = floor(SR*Tf);                           % number of time steps

%% calculate ratio of original length from press postion
for n = 1:6
    if position == -1
        fret = 0;
    else
        fret =1/2^(position/12);
    end
end

%% randomness to improve natureness

famp = famp*((rand(1)+0.1)*0.5+0.4);
dur = dur * round(rand(1)*10+1)/5;

%% fret = 0 mean no note on this string then output a zeros vector
if fret ~= 0 

    %% parameters
    plot_on = 0;                % in-loop plotting on (1) or off (0)

    L = 0.65;                                      % original length (m)
    L_fret = L * fret;                             % actual length (m)

    %% string material
    switch r
        case '013'
            r = [0.00033/2;0.00043/2;0.00066/2;0.00091/2;0.00117/2;0.00142/2];
        case '012'
            r = [0.0003/2;0.00041/2;0.00061/2;0.00081/2;0.00107/2;0.00135/2];
        case '011'
            r = [0.00028/2;0.00038/2;0.00061/2;0.00081/2;0.00107/2;0.00132/2];
    end

    switch material
        case 'steel'
            E = [1.88e11;1.88e11;1.88e11;0.64e11;0.64e11;0.64e11];
            rho = [7690;7950;8220;6930;6610;6540];
        case 'Nylon'
            E = [1.4e9;1.4e9;1.4e9;1.4e9;1.4e9;1.4e9];
            rho = [1150;1150;1150;1150;1150;1150];
        case 'gold'
            E = [7.9e10;7.9e10;7.9e10;7.9e10;7.9e10;7.9e10];
            rho = [19290;19290;19290;19290;19290;19290];
    end
    %T = [3668.1;2058.6;1283.7;727.8;408.5;229.3];
    f0 = [329.63;246.94;195;146.83;110;82.41];   % open strings frequency vector
    T60_1 = [10;12;12;14;15;18];                 % T60 based on frequency 1
    T60_2 = [8;7;6;6;7;5;6];                     % T60 based on frequency 2

    %% physical string parameters

    f1 = 100;                                    % loss based frequency 1
    f2 = 2000;                                   % loss based frequency 2
    omega1 = 2*pi*f1;                            % loss based angular frequency 1
    omega2 = 2*pi*f2;                            % loss based angular frequency 2

    %% I/O
    
    r = r(string_n);                             % radius
    E = E(string_n);                             % young's modular
    rho = rho(string_n);                         % density
    f0 = f0(string_n);                           % fundamental frequency 
    gamma = f0 * 2 * L;                          % calcualte gamma from f0
    T60_1 = T60_1(string_n);                     % T60 of f1
    T60_2 = T60_2(string_n);                     % T60 of f2
    xi = xi * L_fret;                            % coordinate of input
    xo = xo * L_fret;                            % coordinate of output

    %% derived parameters

    A = pi*r^2;                 % string cross-sectional area
    I = 0.25*pi*r^4;            % string moment of inertia
    %T = gamma^2*rho*A;          % Tension for each string
    %c = sqrt(T/(rho*A));        % wave speed;
    %K = sqrt(E*I/(rho*A));      % stiffness constant 
    K = sqrt(B)*(gamma/pi);    % inharmonity effect

    %% loss parameter
    zeta1 = (-gamma^2+sqrt(gamma^4+4*K^2*omega1^2))/(2*K^2);  
    zeta2 = (-gamma^2+sqrt(gamma^4+4*K^2*omega2^2))/(2*K^2);
    sig0 = 6*log(10)*(zeta2/T60_1-zeta1/T60_2)/(zeta2-zeta1);  
    sig1 = 6*log(10)*(-1/T60_1+1/T60_2)/(zeta2-zeta1);         

    %% grid

    k = 1/SR;                                                    % time step
    hmin = sqrt(1/2*(gamma^2*k^2+sqrt(gamma^4*k^4+16*K^2*k^2))); % space step
    N = floor(L_fret/(hmin));                                   % number of space step
    h = L_fret/N;                                                % calculate step again

    %% stability assert
    assert(N<10000,'error! Please check your input parameter'); 
    assert(h >= hmin,'error! Please check your input parameter');
    assert((xi*L)>h),'error! Please check your input parameter';
    assert((L-xi*L)>h,'error! Please check your input parameter');
    assert(dur<Tf,'error! Please check your input parameter');
    assert(r>0||E>0||rho>0||Tf>0,'error! Please check your input parameter');

    %% create force signal

    t_vec = (0:Nf-1)/SR;                    % time vector
    f = zeros(1,Nf);                        % initialise f

    for n = 1:Nf
        if t_vec(n)>=exc_st && t_vec(n)<=exc_st+dur% if in this time window then do
            f(n) = 1/2*famp*(1-cos(itype*pi()*(t_vec(n)-exc_st)/dur)); % signal expression
        end
    end

    %% matrices calulate

    I_mtx = speye(N-1);                                    % eye matrix
    e = ones(N-1,1);                                       % unit vector
    Dxx = 1/h^2*spdiags([e -2*e e], -1:1, N-1,N-1);        % Dxx matrix
    spdiags_1 = spdiags([e -4*e 6*e -4*e e], -2:2, N-1,N-1);
    spdiags_1(1,1) = 5;
    spdiags_1(N-1,N-1) = 5;
    Dxxxx = (1/h^4)*spdiags_1;                             % Dxxxx sparse matrix

    B = 1/(1+sig0*k)*(2 * I_mtx + gamma^2*k^2*Dxx -K^2*k^2*Dxxxx+2*sig1*k*Dxx);% coefficient of u1
    C = ((-1+sig0*k)*I_mtx-2*sig1*k*Dxx)/(1 + sig0*k);                         % coefficient of u2

    %% calcualte input/ouput coordinate index
    li = floor(xi/h);                                      % index of input position
    lo = floor(xo/h);                                      % index of output position

    J = zeros(N-1,1);                                      % initialse J
    J(li) = 1;
    J = k^2/(rho*A*h)*J;                                   % calculate J
    
    c_v = zeros(N-1,1);
    c_v(lo) = 1;                                           % output selection

    %% initialise scheme variables

    u2 = zeros(N-1,1);              % state
    u1 = u2;                        % state
    u = u2;                         % state

    y = zeros(Nf,1);                % output

    %% main loop

    tic
    for n=1:Nf

        % update state, and insert current value of f.

        u = B*u1+C*u2+J*f(n);
        
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
else
    y = zeros(Nf,1);
end
%% play sound

%soundsc(y,SR);

%% plot spectrum

figure(2)

yfft = 10*log10(abs(fft(y)));
plot([0:Nf-1]'*(SR/Nf), yfft, 'k')
xlim([0 SR/2])
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Transform of output')
end
