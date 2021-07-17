clc;
clear all;
close all;
SR = 44100;
xi = 0.8;                                 % coordinate ratio of input
xo = 0.7;                                 % coordinate ratio of output
famp = 1;                                 % input amplitude
dur = 0.0001;                             % input duration
exc_st = 0.1;                             % input timing
itype = 1;                                % stike type
speed = 0.1;                              % speed to play upward or downward
B = 0.0001;                               % inharmonity parameter
Tf = 8;                                   % output duration
r = '013';                                % select global string stickness from 011/012/013
material = 'Nylon';                       % select material from steel/Nylon/gold

%% chord

y = BTB_s2119032_Lai_freqdependent_chord('Cmaj7','down',xi,xo,famp,dur,exc_st,itype,'011','steel',B,speed);
soundsc(y,SR);

% y = BTB_s2119032_Lai_freqdependent_chord('Bbm7','up',xi,xo,famp,dur,exc_st,itype,r,material,B,0.01);
% soundsc(y,SR);

%% open strings panning

% y  =  BTB_s2119032_Lai_freqdependent_note_sequency(xi,xo,famp,dur,exc_st,itype,r,material,B,0,Tf);
% soundsc(y,SR);

%% single note
% y =  BTB_s2119032_Lai_freqdependent(xi,xo,famp,dur,exc_st,6,itype,'013','Nylon',B,0,Tf);
% soundsc(y,SR);
