%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:BTB_s2119032_Lai_freqdependent_chord
%%%         Author:Nuolin Lai
%%%         Create Date:11/03/2021
%%%         Last modify date:11/03/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y  =  BTB_s2119032_Lai_freqdependent_chord(chord,direct,xi,xo,famp,dur,exc_st,itype,r,material,B,speed)

%% chord

switch chord
    case 'A7'
        seq = [-1;13;12;11;12;-1];
    case 'Am7'
        seq =[-1;5;5;5;-1;5];
    case 'Bb7'
        seq = [-1;8;7;6;-1;6];
    case 'Bbm7'
        seq = [-1;2;1;3;1;-1];
    case 'B7'
        seq = [-1;9;8;7;-1;7];
    case 'Bm7'
        seq=[-1;6;7;7;-1;7];
    case 'C7'
        seq = [-1;10;9;8;-1;8];
    case 'Cm7'
        seq = [-1;4;3;5;3;-1];
    case 'Cmaj7'
        seq = [-1;5;4;5;3;-1];
    case 'Db7'
        seq = [-1;11;10;9;-1;9];
    case 'Dbm7'
        seq = [-1;6;5;6;4;-1];
    case 'D7'
        seq = [-1;7;5;7;5;-1];
    case 'Dmaj7'
        seq = [-1;7;5;7;5;-1];
    case 'Dm7'
        seq = [-1;6;5;7;5;-1];
    case 'Eb7'
        seq = [-1;4;6;5;6;-1];
    case 'Ebm7'
        seq = [-1;8;7;8;6;-1];
    case 'E7'
        seq = [-1;8;7;6;7;-1];
    case 'F7'
        seq = [-1;1;2;1;-1;1];
    case 'Fm7'
        seq = [-1;9;8;10;8;-1];
    case 'F#7'
        seq = [-1;9;9;8;9;-1];
    case 'G7'
        seq = [-1;3;4;5;-1;3];
    case 'Gmaj7'
        seq = [-1;3;4;4;-1;3];
    case 'G#7'
        seq = [-1;11;11;10;11;-1];
end


%% parameters
SR= 44100;
Nf = floor(SR*5);
interval_N = floor(speed*SR);
y = zeros(Nf+6*interval_N,1);

for n = 1:6
    switch direct
        case 'down'
            index_n = 6-n;
        case 'up'
            index_n = n;
    end

    y1 = BTB_s2119032_Lai_freqdependent(xi,xo,famp,dur,exc_st,n,itype,r,material,B,seq(n),5);
    y(1+index_n*interval_N:Nf+index_n*interval_N,1) = y(1+index_n*interval_N:Nf+index_n*interval_N,1) + y1;

end

