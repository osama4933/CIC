function [Data_buff peak_amp Up_ind FFO] = dnsamp_buff(Data_stack,Upchirp_ind)
%dnsamp_buff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data_stack contains All possible downsampled Buffers with different starting sample for downsampling
% e.g.  Signal = [ 1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16]
% ^ this Signal can be Downsampled by a factor of 8(upsampling_factor) in
% following ways
% [1    9], [2      10], [3     11], [4     12], [5     13], [6     14], [7     15], [8     16]
% Now out of 8 data_stack possibilities, there exists one buffer that gives
% very good quality frequency Tracks (FFT followed by Dechirping gives the most crispy Peak)
% Definition of Good Quality Frequency Track: Energy of FFT peak does not
% leaks into adjacent bins.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load Parameters
SF = param_configs(1);
BW = param_configs(2);
Fs = param_configs(3);
N = 2^SF;
num_preamble = param_configs(4);
num_sync = param_configs(5);
num_DC = param_configs(6);

DC = conj(sym_to_data_ang([1],N));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Compute and Correct Frequency Offsets for each Preamble Detected in each Data_stack and Find the Peak Statistics needed for demodulation
Up_ind = [];
peak_amp = [];
Data_buff = [];
FFO = [];
ffo = [];
% n_pnt is the fft Factor - fft(Signal, n_pnt * length(Signal))
n_pnt = 16;             


% iterate over all Upchirps that qualified 8 consecutive Peak condition
for k = 1:size(Upchirp_ind,1)
    if(Upchirp_ind(k,1) - N <= 0)
        continue;
    end
    close all
    in = [];
    % iterate overall downsampled buffers
    for m = 1:size(Data_stack,1)
            data_wind = [];
            data_fft = [];
            freq_off = [];
            % ind_temp contains the Frequency Bins around bin 1 where a
            % Preamble Peak can lie
            ind_temp = [1:5*n_pnt (N*n_pnt)-(4*n_pnt):(N*n_pnt)];
            % iterate over all Preambles
            for j = 1:num_preamble
                data_wind = Data_stack(m,Upchirp_ind(k,1) : Upchirp_ind(k,1) + (num_preamble*N) -1);
                data_fft(j,:) = abs(fft(data_wind((j-1)*N + 1:j*N) .* DC(1:N),n_pnt*N));
                [~,c(j)] = max(data_fft(j,ind_temp));
                c(j) = ind_temp(c(j));
                % Handle -ve and +ve Frequency Offsets Accordingly
                if(c(j) > (n_pnt*N)/2)
                    freq_off = [freq_off ( (N*n_pnt) - c(j) ) / n_pnt];     % +ve offset
                else
                    freq_off = [freq_off -1*( c(j) - 1 ) / n_pnt];          % -ve offset
                end
            end
            % average the frequency offset of 6 middle Preambles
            freq_off = sum( freq_off(2:num_preamble-1) ) / (num_preamble - 2);
            ffo = [ffo freq_off];
            % Correct for the Frequency Offset in corresponding Data_Stack
            Data_freq_off(m,:) = Data_stack(m,:) .* exp( (1i*2*pi*(freq_off./N)) .* (1:length(Data_stack(m,:))) );
            
            clear data_wind data_fft ind_temp
            % ind_temp contains the Frequency Bins around bin 1 where a
            % Preamble Peak can lie, assumption (-5*BW/2^SF <= Freq_off <= 5*BW/2^SF)
            ind_temp = [1:5 (N-4):N];
            a = [];
            % for the frequency offset corrected Data Stack, find FFT of Preamble to get Peak Statistics 
            for j = 1:num_preamble
                data_wind = Data_freq_off(m,Upchirp_ind(k,1) : Upchirp_ind(k,1) + (num_preamble*N) -1);
                data_fft(j,:) = abs(fft(data_wind((j-1)*N + 1:j*N) .* DC(1:N),N));
                [a(j),c(j)] = max(data_fft(j,ind_temp));
                c(j) = ind_temp(c(j));
            end
            peak_stats(k,m,1) = mean(a);
            peak_stats(k,m,2) = var(a);
            peak_stats(k,m,3) = std(a);
            
            %%  Find the Right Data_stack to work with
            % first find the stft of given stack at the Preamble Region,
            % Spec is a 2D Matrix, rows - Freq. Bins & Col. - Time Samples
            Spec = stft_v1(Data_freq_off(m,Upchirp_ind(k,1) - N:Upchirp_ind(k,end) + N - 1 - N),N,DC(1:N),0,0);
            temp = [];
            freq_track_qual = [];
            pream_peak_ind = [];
            adj_ind = [];
            % row_ind contains the Frequency Rows around bin 1 where a
            % Preamble Peak can lie
            row_ind = [N-5:N 1:6];
            count = 1;
            for i = row_ind
                temp(count) = sum(abs(Spec(i,:)));
                count = count + 1;
            end
            % Frequency Track in row containing Preamble should have
            % maximum energy
            [~,ind] = max(temp);
            pream_peak_ind = row_ind(ind);
            % Find row indices for Preamble row + 1 & - 1
            adj_ind = [mod(pream_peak_ind-1,N) mod(pream_peak_ind+1,N)];
            if(sum(adj_ind == 0) == 1)
                adj_ind(find(adj_ind == 0)) = N;
            end
            % A good quality frequency track for a preamble is one that has
            % least energy leakage in adjacent rows (this promises very sharp FFT peaks)
            freq_track_qual = ( sum(abs(Spec(pream_peak_ind,:))) - sum(abs(Spec(adj_ind(1),:))) ) + ( sum(abs(Spec(pream_peak_ind,:))) - sum(abs(Spec(adj_ind(2),:))) );
            in = [in freq_track_qual];
    end
    % choosing the best Data_stack based on maximum energy difference from
    % adjacent bins
    [~,b] = max(in);
    % output frequency offset corrected buffer with relevant, Peak
    % statistics and frequency offsets
    Data_buff = [Data_buff; Data_freq_off(b,:)];
    FFO = [FFO; ffo(b)];
    peak_amp = [peak_amp; reshape(peak_stats(k,b,:),1,[])];
    Up_ind = [Up_ind; Upchirp_ind(k,:)];
end

end
