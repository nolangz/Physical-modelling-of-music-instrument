%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         PBMMI:Exercise2
%%%         Author:Nuolin Lai
%%%         Create Date:18/01/2020
%%%         Last modify date:18/01/2020
%%%         General High order FIR filter 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y] = myconv(h,x)


Nb = length(h);
Nx = length(x);

%change x to N-length Kronecker delta function 
%x = [1;zeros(Nx-1,1)]

%create impulse response vector h
%Ny = Nx+Nb;

%frequency vector
% fx = (0:Nx-1)';
% fy = (0:Ny-1)';
Y = 1;

%fft
X = fft(x,max(Nb,Nx));
H = fft(h,max(Nb,Nx));

Y = X.*H;

y= ifft(Y,'symmetric');


%xpad(Nb+1-m:end-m)
% %DFT
% XF = fft(x);
% YF = fft(y);
% L1 = length(XF);
% plot(fy,20*log10(abs(YF)));
% 
% T = fvtool(h,1);
% 
% %play the sound
% sound(y,Fs)
