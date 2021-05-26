function [Preamble_ind, bin_offsets, Data_out, Peak_amp, ffo] = filter_false_postives(Data_stack,Upchirp_ind,Peak,FFO)
%FILTER_FALSE_POSTIVES Summary of this function goes here
% This function takes in detected Possible Preambles and there relevant
% data_stacks and then looks at points of concern for the presence of 2
% SYNC-WORDS. This final filter removes false Positives with a very high
% certanity

% load parameters
SF = param_configs(1);
BW = param_configs(2);
Fs = param_configs(3);
N = 2^SF;
num_preamble = param_configs(4);
num_sync = param_configs(5);
num_DC = param_configs(6);
S1 = param_configs(12);
S2 = param_configs(13);

DC = conj(sym_to_data_ang(ones(1,num_preamble),N));	% Generate a num_preamble # of downchirps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ffo = [];
Preamble_ind = [];
bin_offsets = [];
Data_out = [];
Peak_amp = [];

% row_ind contains the Frequency Bins around bin 1 where a
% Preamble Peak can lie, assumption (-6*BW/2^SF <= Freq_off <= 6*BW/2^SF)
row_ind = [N - 4:N + 1 1:6];
for m = 1:size(Upchirp_ind,1)

    % extract 8 Preamble long window
    data_wind = Data_stack(m,Upchirp_ind(m,1) : Upchirp_ind(m,1) + ((num_preamble) * N) - 1);
    close all

    % Compute STFT to accurately find the Preamble Frequency Bin and any
    % bin offset (if any left)
    [Spec] = stft_v2(data_wind.*DC,N);
    temp = [];
    count = 1;
    for i = row_ind
        temp(count) = sum(abs(Spec(i,:)));
        count = count + 1;
    end
    [~,ind] = max(temp);
    pream_peak_ind(m) = row_ind(ind);

    % from the Preamble find expected SYNC WORD bins
    sync1_ind = mod(pream_peak_ind(m) + S1,N);
    sync2_ind = mod(pream_peak_ind(m) + S2,N);
    if(sync1_ind == 0)
        sync1_ind = N;
    end
    if(sync2_ind == 0)
        sync2_ind = N;
    end
    
    % Extract windows corresponding to 2 SYNC-WORDS
    sync_wind = Data_stack(m,Upchirp_ind(m,num_preamble) + N : Upchirp_ind(m,num_preamble) + N + (num_sync*N) - 1);
    % compute the thresholds for SYNC WORD's Peak and perform dechirping + FFT
    sync_threshold_up = Peak(m,1) + 0.5 * Peak(m,1);
    sync_threshold_low = Peak(m,1) - 0.5 * Peak(m,1);
    
    sync_word1 = abs(fft(sync_wind(1:N).*DC(1:N)));
    sync_word2 = abs(fft(sync_wind(N+1:end).*DC(1:N)));
    
    if(sync_threshold_low < (2 * sum(sync_word1) / N))
        sync_threshold_low = (2 * sum(sync_word1) / N);
    elseif( sync_threshold_low < (2 * sum(sync_word2) / N))
        sync_threshold_low = (2 * sum(sync_word2) / N);
    end
    
    % Extract Peaks qualifying Peak Thresholds
    syn1_pnts = get_bounded_max(sync_word1,sync_threshold_up,sync_threshold_low);
    syn2_pnts = get_bounded_max(sync_word2,sync_threshold_up,sync_threshold_low);

    % Check If both SYNC WORDS are present
    if(sum(syn1_pnts == sync1_ind) && sum(syn2_pnts == sync2_ind))
        Preamble_ind = [Preamble_ind; Upchirp_ind(m,:)];
        if(pream_peak_ind(m) < N / 2)
            bin_offsets = [bin_offsets 1 + (-mod(pream_peak_ind(m),N))];
        else
            bin_offsets = [bin_offsets mod(N + 2 - pream_peak_ind(m),N)];
        end
        Data_out = [Data_out; Data_stack(m,:)];
        Peak_amp = [Peak_amp; Peak(m,:)];
        ffo = [ffo; FFO(m)];
    end
end


end

