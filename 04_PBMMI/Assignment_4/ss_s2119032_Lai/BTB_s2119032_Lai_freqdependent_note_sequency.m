%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:BTB_s2119032_Lai_freqdependent_note_sequency
%%%         Author:Nuolin Lai
%%%         Create Date:13/03/2021
%%%         Last modify date:13/03/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y  =  BTB_s2119032_Lai_freqdependent_note_sequency(xi,xo,famp,dur,exc_st,itype,r,material,B,position,Tf)

%% parameters
SR= 44100;
Nf = floor(SR*Tf);
speed = 0.5;
interval_N = floor(speed*SR);
y = zeros(Nf+6*interval_N,2);

%panning
pan = [-1;0.8;-0.6;0.4;-0.2;0];
pan = flip(pan);
panL = -0.5*pan+0.5;
panR = 0.5*pan+0.5;

%sum
for n = flip(1:6)
    y1 = BTB_s2119032_Lai_freqdependent(xi,xo,famp,dur,exc_st,n,itype,r,material,B,position,Tf);
    y(1+(6-n)*interval_N:Nf+(6-n)*interval_N,1) = y(1+(6-n)*interval_N:Nf+(6-n)*interval_N,1) +panL(n)*y1;
    y(1+(6-n)*interval_N:Nf+(6-n)*interval_N,2) = y(1+(6-n)*interval_N:Nf+(6-n)*interval_N,2) +panR(n)*y1;
end
end
