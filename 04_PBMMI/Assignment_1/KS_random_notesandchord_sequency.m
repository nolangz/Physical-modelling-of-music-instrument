%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:KS_random_notesandchord_sequency
%%%         Author:Nuolin Lai
%%%         Create Date:30/01/2021
%%%         Last modify date:30/01/2021
%%%         Random music based on C D E F G A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;

%% global parameter
fA = 440;                           %standard frequency in Hz
fC = fA*2^(3/12)/4;                 %frequency of low C in guitar
syn_min = 1/2;                      %Minimal syncopated note
pitch_c = [0 4 5 7 9];              %pitch of chord
pitch_n = [0 2 4 5 7 9];            %pitch of notes
A = [1, 2, 4];                      %syncopated notes

%% input
BPM = 70;                           %beats per minute
Fs = 44100;                         %sampling rate
NBar = 8;                           %number of bar
Mu= 1/2;                            %position of pick
R = 0.95;                           %dynamics parameter R
S = 0.5;                            %the decay stretching factor

%% calculate parameter
Nb = 60/BPM*Fs;                     %number per beats
fre_n   = fC*2.^(pitch_n/12);       %frequency of notes
random_c = pitch_c(randi(numel(pitch_c),8,1));%chord sequency
random_n = pitch_n(randi(numel(pitch_n),2*NBar*4,1));%chord sequency
A = syn_min*A;                      %actual syncopated notes
Nn = length(random_n);              %initialise number of notes
N0 = zeros(Nn,1);                   %initialise sample number of notes

%% create actual notes based on chords length
N = 0;                              %initialise N
for p = 1:Nn
    random_num = A(randi(numel(A),1,1));
    N0 (p) = random_num*Nb;
    N = N + N0(p);
    if N >= Nb*NBar*4
        break
    end
end
N0 = N0(N0~=0);
N0 = floor(N0);
Nn_exact = length(N0);

%% chords bank
y0 = KS_chord_acoustics(BPM,4,pitch_c(1),0,0,Fs,R,S,Mu);
L = length(y0);
y = zeros(L,10);
y(:,0+1)=KS_chord_acoustics(BPM,4,0,0,0,Fs,R,S,Mu);
y(:,4+1)=KS_chord_acoustics(BPM,4,4,1,0,Fs,R,S,Mu);
y(:,5+1)=KS_chord_acoustics(BPM,4,5,0,0,Fs,R,S,Mu);
y(:,7+1)=KS_chord_acoustics(BPM,4,7,0,0,Fs,R,S,Mu);
y(:,9+1)=KS_chord_acoustics(BPM,4,9,1,0,Fs,R,S,Mu);

%% notes bank
x0 = KS_singlenotes_acoustics(BPM,4,fre_n(1),Fs,R,S,Mu);
M  = length(x0);
x  = zeros(M,10);
x(:,0+1) = KS_singlenotes_acoustics(BPM,4,fre_n(1),Fs,R,S,Mu);
x(:,2+1) = KS_singlenotes_acoustics(BPM,4,fre_n(2),Fs,R,S,Mu);
x(:,4+1) = KS_singlenotes_acoustics(BPM,4,fre_n(3),Fs,R,S,Mu);
x(:,5+1) = KS_singlenotes_acoustics(BPM,4,fre_n(4),Fs,R,S,Mu);
x(:,7+1) = KS_singlenotes_acoustics(BPM,4,fre_n(5),Fs,R,S,Mu);
x(:,9+1) = KS_singlenotes_acoustics(BPM,4,fre_n(6),Fs,R,S,Mu);

%% create chord sequency
y2 = zeros(Nb*NBar*4+L,1);
for i = 1:NBar
    for j = 1:4   
        y2(((i-1)*4+j-1)*Nb+1:((i-1)*4+j-1)*Nb+L,1)=y2(((i-1)*4+j-1)*Nb+1:((i-1)*4+j-1)*Nb+L,1)+y(:,random_c(i)+1);
    end
end

%% create notes
x2 = zeros(sum(N0)+M,1);
for q = 1:Nn_exact
    if q == 1
        x2(1:M,1)= x2(1:M,1)+x(:,random_n(q)+1);
    else
    x2(sum(N0(1:q))+1:sum(N0(1:q))+M,1)=x2(sum(N0(1:q))+1:sum(N0(1:q))+M,1)+x(:,random_n(q)+1);
    end
end

%% Normalise
s = x2(1:Nb*NBar*4+L)+y2*0.75;
s = s/max(abs(s));
sound(s,Fs)

        
    

