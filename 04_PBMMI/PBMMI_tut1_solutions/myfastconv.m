% -------------------------------------------------
% FFT COVOLUTION
% Alistair Carson 2021
% -------------------------------------------------
function y = myfastconv(h,x)
   Nx = length(x);
   Nb = length(h);
   Nfft = 2^(nextpow2(Nx + Nb));
   y=ifft(fft(h,Nfft).*fft(x,Nfft),'symmetric');
end

