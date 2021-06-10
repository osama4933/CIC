function [Spec] = stft_v2(Buffer,f_n)
% This function produces a spectrogram of Dechirped LoRa signal and if plotted, gives a visualization of
% continuous frequency tracks. (Best Frequerncy Resolution, worst Time Resolution (we already have that from Preamble detection))

w_n = f_n;
Buff_len = length(Buffer); % Signal length
% Frequency axis
f_n = ceil(f_n/2) * 2+1;
Lf = (f_n - 1)/2;
% Time axis
w_n = ceil(w_n/2) * 2+1;
Lw = (w_n - 1)/2;
% Initialize Spectrum to zero with appropriate size
Spec = zeros(f_n,Buff_len); 
% Sliding window over signal
for iter = 1:Buff_len
    i_l = min([iter-1, Lw, Lf]);
    i_r = min([Buff_len-iter, Lw, Lf]); 
    iter_ind = -i_l:i_r;
    ind1 = iter_ind + iter;   % Time Indexing of the original signal
    ind = iter_ind + Lf +1;     % Frequency Indexing of the martix
    
    temp_buff = Buffer(ind1);
    Spec(ind, iter) = temp_buff;
end
% Computing FFT of stacked Windows
Spec = fft(Spec);
Spec = (Spec*2) / f_n;  % normalizing the FFTs
end

