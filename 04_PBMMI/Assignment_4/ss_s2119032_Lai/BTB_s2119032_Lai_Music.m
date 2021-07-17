%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:Stiffness string - random music
%%%         Author:Nuolin Lai
%%%         Create Date:13/03/2021
%%%         Last modify date:13/03/2021
%%%         include:
%%%         1. random music based on parameter and sequence of note/chord
%%%         2. random panning for each notes
%%%         3. random upward or downward for playing chord
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;

%% global bank

updown = ["up","down"];                   % play upward or downward
chord = ["A7","Am7","Bb7","Bbm7","B7","Bm7","C7","Cm7","Cmaj7","Db7","Dbm7","D7","Dmaj7","Dm7","Eb7","Ebm7","E7","F7","Fm7","F#7","G7","Gmaj7","G#7"];

%% input parameter to create random 

type_c = [1 5 10 13 17 20 23];            % index of chord
note_string = [2,3,2,4,3,4,3,1,1,2,1];    % select a sequence of notes - string number
note_fret_n = [9,6,7,9,9,7,7,9,7,10,10];  % select a sequence of notes - fret number
note_index = 1:length(note_string);       % index of fret number
chord_index = 1:length(type_c);           % index of chord
A = [1, 2, 4];                            % syncopated notes
syn_min = 1/4;                            % Minimal syncopated note
BPM = 120;                                % beats per minute
SR = 44100;                               % sampling rate
Nb = 60/BPM*SR;                           % number per beats
NBar = 8;                                 % number of bar
C_inbar = 4;                              % chord number in one bar      
SR = 44100;                               % sample rate
xi = 0.8;                                 % coordinate ratio of input
xo = 0.7;                                 % coordinate ratio of output
famp = 1;                                 % input amplitude
dur = 0.0001;                             % input duration
exc_st = 0.1;                             % input timing
itype = 1;                                % stike type
speed = 0.01;                             % speed to play upward or downward
B = 0.0001;                               % inharmonity parameter
r = '011';                                % select global string stickness from 011/012/013
material = 'steel';                       % select material from steel/Nylon/gold
%% calculate parameter
pan = (rand(2*NBar*C_inbar,1)-0.5)*2;     % random pan for each note
panL = -0.5*pan+0.5;                      % panning L channel volume for each note
panR = 0.5*pan+0.5;                       % panning R channel volume for each note
random_c = chord_index(randi(numel(chord_index),8,1));           % chord sequency index
random_n = note_index(randi(numel(note_index),2*NBar*C_inbar,1));% note sequency index
random_ud = round(rand(1,8)*1)+1;                                % random upward or downward
A = syn_min*A;                            % actual syncopated notes
Nn = length(random_n);                    % initialise number of notes
N0 = zeros(Nn,1);                         % initialise sample number of notes

%% create actual notes based on chords length

% to let chords and melody end at the almost same time

N = 0;                              %initialise N
for p = 1:Nn
    random_num = A(randi(numel(A),1,1));
    N0 (p) = random_num*Nb;
    N = N + N0(p);
    if N >= Nb*NBar*C_inbar
        break
    end
end
N0 = N0(N0~=0);
N0 = floor(N0);
Nn_exact = length(N0);              % output the actual notes number 

%% chords bank
% initialse chords
y0 = BTB_s2119032_Lai_freqdependent_chord(chord(type_c(1)),'down',xi,xo,famp,dur,exc_st,itype,r,material,B,speed);
L = length(y0);
y = zeros(L,7);
% initialse chord bank
for n = 1:length(type_c)
    y(:,n) = BTB_s2119032_Lai_freqdependent_chord(chord(type_c(n)),updown(random_ud(n)),xi,xo,famp,dur,exc_st,itype,r,material,B,speed);
end

%% notes bank
% initialse notes
x0 = BTB_s2119032_Lai_freqdependent(xi,xo,famp,dur,exc_st,note_string(1),itype,r,material,0.0001,note_fret_n(1),5);
M  = length(x0);
x  = zeros(M,10);
% initialse notes bank
for n = 1:11
    x(:,n) = BTB_s2119032_Lai_freqdependent(xi,xo,famp,dur,exc_st,note_string(n),itype,r,material,0.0001,note_fret_n(n),5);
end


% create chord sequency output
y2 = zeros(Nb*NBar*C_inbar+L,1);        % initialise output vector
empty = [1 2 3 4];                      % random bank
for i = 1:NBar 
    e = empty(randi(numel(empty),1,1)); % random one beats to be silence in each bar
    for j = 1:C_inbar  
        if j ~= e                       % when beats number equal to silence beats number, add a zeros vector into output
            y2(((i-1)*C_inbar+j-1)*Nb+1:((i-1)*C_inbar+j-1)*Nb+L,1)=y2(((i-1)*C_inbar+j-1)*Nb+1:((i-1)*C_inbar+j-1)*Nb+L,1)+y(:,random_c(i));
        else
            y2(((i-1)*C_inbar+j-1)*Nb+1:((i-1)*C_inbar+j-1)*Nb+L,1)=y2(((i-1)*C_inbar+j-1)*Nb+1:((i-1)*C_inbar+j-1)*Nb+L,1)+zeros(L,1);
        end
    end
end

%% create notes
x2 = zeros(sum(N0)+M,2);
for q = 1:Nn_exact
    if q == 1 % put the first note into output vector
        x2(1:M,1)= x2(1:M,1)+x(:,random_n(q));
        x2(1:M,2)= x2(1:M,2)+x(:,random_n(q));
    else      % put the rest of notes into output vector
    x2(sum(N0(1:q))+1:sum(N0(1:q))+M,1)=x2(sum(N0(1:q))+1:sum(N0(1:q))+M,1)+panL(q) * x(:,random_n(q));
    x2(sum(N0(1:q))+1:sum(N0(1:q))+M,2)=x2(sum(N0(1:q))+1:sum(N0(1:q))+M,2)+panR(q) * x(:,random_n(q));
    end
end

%% Normalise
% compare the sample nubmer of chords output and notes output and combine them
y2(:,2) = y2(:,1);
if length(y2(:,1)) < length(x2(:,1))  
    s(:,1) = x2(1:length(y2),1)+y2(:,1)*0.75;
    s(:,2) = x2(1:length(y2),2)+y2(:,2)*0.75;
else
    s(:,1) = x2(:,1) + y2(1:length(x2),1)*0.75;
    s(:,2) = x2(:,2) + y2(1:length(x2),2)*0.75;
end

%% play
s = s/max(max(abs(s)));
soundsc(s,SR)


        
    

