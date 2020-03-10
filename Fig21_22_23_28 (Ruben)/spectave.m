function [f,spm] = spectave(signal,fs,NFFT)
%Compute average spectrum of a signal(:,i) over different runs i at a
%sampling frequency fs
[a,b] = size(signal);
for i=1:b
    [sp(:,i),f] = compspe2(detrend(signal(:,i)),fs,NFFT);
end
spm = mean(sp,2);
end

