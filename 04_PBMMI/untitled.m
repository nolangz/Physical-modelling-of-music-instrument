T = 5;
Fs = 44100;
N = T*Fs;
t = (0:N-1)/Fs;
f1 = 600;
f2 = 5;

o2 = 0.35*cos(ppa(0.35*cos(2*pi*f1*t)+2*pi*f2*t));
phase = 2*pi*(sin(2*pi*f1*t)+2*pi*f2*t);
%% plot waveform and spectrum
subplot(2,1,1);
plot(t,o2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Waveform')
subplot(2,1,2);
spe = 20*log10(abs(fft(o2)));
freq = (0:N-1)/N*Fs;
semilogx(freq,spe);
xlabel('frequency (Hz)');
ylabel('Amplitude(dB)');
title('Spectrum')
ylim([-20 100]);
xlim([20 Fs/2]);
sound(o2,Fs);
