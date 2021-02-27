%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:KS_chord_acoustics
%%%         Author:Nuolin Lai
%%%         Create Date:28/01/2021
%%%         Last modify date:30/01/2021
%%%         Input BPM,syncopated note type(1,1/2,1/4),fundamental freqeuncy
%%%         Sampling frequency,Dynamics parameter,decay stretching factor
%%%         BPM = 60~120 (suggestion)
%%%         syn = 4,2,1,1/2,1/4,1/8....
%%%         chord = 0 1  2 3  4 5 6  7 8  9 10 11 corresponse to 
%%%                 C C# D D# E F F# G G# A Bb   
%%%         type: 'Major'=0 or 'Minor'=1
%%%         octave:0 or 1
%%%         Fs: sampling rate
%%%         R（0，1）: dynamics parameter R
%%%         S（0，1）: the decay stretching factor
%%%         Mu (0,1): pick position (0,1) 1/2 is in the middle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function y = KS_chord_acoustics(BPM,syn,chord,type,octave,Fs,R,S,Mu);
%create chord notes frequency vector
fA = 440;                           % standard frequency in Hz
fC = fA*2^(3/12)/4*2^(octave);     % frequency of low C in guitar
fc = fC/2*2^(octave);              % C notes out of guitar range
if chord >= 4 && chord <= 8
    if type == 0
        p0 = chord;                 %6th string
        p1 = p0+7;                  %5th string
        p2 = p1+5;                  %4th string
        p3 = p2+4;                  %3rd string
        p4 = p3+3;                  %2nd string
        p5 = p4+5;                  %1st string
        p  = [p0 p1 p2 p3 p4 p5];   %pitch vector
        f  = fc*2.^(p/12);          %frequency vector
    elseif type == 1
        p0 = chord;                 %6th string
        p1 = p0+7;                  %5th string
        p2 = p1+5;                  %4th string
        p3 = p2+2;                  %3rd string
        p4 = p3+5;                  %2nd string
        p5 = p4+5;                  %1st string
        p  = [p0 p1 p2 p3 p4 p5];   %pitch vector
        f  = fc*2.^(p/12);          %frequency vector 
    end
elseif chord >= 0 && chord <= 3
    if type == 0
        p0 = chord;                 %5th string
        p1 = p0+7;                  %4th string
        p2 = p1+5;                  %3rd string
        p3 = p2+4;                  %2nd string
        p4 = p3+3;                  %1st string
        p  = [p0 p1 p2 p3 p4];      %pitch vector
        f  = fC*2.^(p/12);          %frequency vector
    elseif type == 1
        p0 = chord;                 %5th string
        p1 = p0+7;                  %4th string
        p2 = p1+5;                  %3rd string
        p3 = p2+3;                  %2nd string
        p4 = p3+4;                  %1st string
        p  = [p0 p1 p2 p3 p4];      %pitch vector
        f  = fC*2.^(p/12);          %frequency vector  
    end
else
    if type == 0
        p0 = chord;                 %5th string
        p1 = p0+7;                  %4th string
        p2 = p1+5;                  %3rd string
        p3 = p2+4;                  %2nd string
        p4 = p3+3;                  %1st string
        p  = [p0 p1 p2 p3 p4];      %pitch vector
        f  = fc*2.^(p/12);          %frequency vector
    elseif type == 1
        p0 = chord;                 %5th string
        p1 = p0+7;                  %4th string
        p2 = p1+5;                  %3rd string
        p3 = p2+3;                  %2nd string
        p4 = p3+4;                  %1st string
        p  = [p0 p1 p2 p3 p4];      %pitch vector
        f  = fc*2.^(p/12);          %frequency vector  
    end   
end
f0 = fc*2^(p0/12);          %fundamental frequency of  6th string.
y0 = KS_singlenotes_acoustics(BPM,syn,f0,Fs,R,S,Mu);
L = length(y0);
y  = zeros(L,1);            %initialise y
for f0 = f
    y = y + KS_singlenotes_acoustics(BPM,syn,f0,Fs,R,S,Mu);
end
y = y/max(y);               %Normalise
end