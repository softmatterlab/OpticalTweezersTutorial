 function [psd, absFFTout, getfft, f]=FFTnormalized(signal, fsample )

NFFT = 2^nextpow2(size(signal,1)); % Next power of 2 from length of signal
getfft = fft(signal,NFFT)/size(signal,1);
absFFTout = 2*abs(getfft(1:NFFT/2+1,:)).^2;  % get the positive frequencies 
psd = absFFTout*size(signal,1)/fsample;
f = fsample/2*linspace(0,1,NFFT/2+1);  % get the frequency axis