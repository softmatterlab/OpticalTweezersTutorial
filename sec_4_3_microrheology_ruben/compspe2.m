function [sp,f]=compspe2(AA,fs,NFFT)
%power spectral density
%fs=sampling frequency
aa = length(AA);
OVERLAP=round(aa*0.1);
%[sp,f]=irf_psd(AA,NFFT,fs,hanning(NFFT),OVERLAP); 
[sp,f]=pwelch(AA,[],OVERLAP,NFFT,fs); 


return
