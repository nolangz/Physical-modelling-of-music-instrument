%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:Exercise1
%%%         Author:Nuolin Lai
%%%         Create Date:18/01/2020
%%%         Last modify date:18/01/2020
%%%         compare conv function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;

%import audio
[x,Fs]   = audioread('birchcanoe.wav');

%create impulse response vector h
h = [1;1;1];
x = x';

% %compare myconv and conv
% y1 = myconv(h,x);
% y2 = conv(h,x);
% comparison = y1 - y2;
% % done

N= 10000;
Nb = N-1;
n=0:Nb;
h = 0.5-0.5*cos(2*pi*n/(Nb+1));


tic;
for t = 1:10
    y1 = myfastconv(h,x);
end
T1 = toc;


tic;
for t = 1:10
    y2 = conv(h,x);
end
T2 = toc;