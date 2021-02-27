%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:KS random_chord_sequency
%%%         Author:Nuolin Lai
%%%         Create Date:30/01/2020
%%%         Last modify date:30/01/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
%% global parameter
a = round(rand(1,30)*11) ;          %random number
syn_min = 1/4;                      %minimal syncopated note
A = [1, 2, 4, 8];                   %syncopated note range
A = syn_min*A;                      %actual syncopated note range


%% input
BPM = 70;                           %beats per minute
Fs  = 44100;                        %sampling rate
Mu  = 1/2;                          %position of pick
R   = 0.95;                         %dynamics parameter R
S   = 0.5;                          %the decay stretching factor
Nchord = 30;                        %number of notes
type_c = 1;                         %Major or Minor chord

%% create chord sequency
N = 0;                              % initialise N
syn = zeros(1,30);                  % initialise syncopated note
for i = 1:Nchord
        random_num = A(randi(numel(A),1,1));
        N0 = random_num*60/BPM*Fs;
        N = N + N0;
        syn(i) = random_num;
        y = KS_chord_acoustics(BPM,syn(i),a(i),1,0,Fs,R,S,Mu);
        sound(y,Fs);
end

        
    

