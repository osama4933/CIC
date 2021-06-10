function [Upchirp_ind] = UC_location_corr_DC_based(Data,DC_ind)
%UC_location_corr_DC_based 
% this function takes in the Data buffer, the Downchirp indices from
% previous function and looks locally for presence of an Upchirp in order
% to increase the certainity of the presence of a pkt

% chirp variables
SF = param_configs(1);
BW = param_configs(2);
Fs = param_configs(3);
N = 2^SF;
upsampling_factor = Fs/BW;

% LORA pkt variables
num_preamble = param_configs(4);
num_sync = param_configs(5);
num_DC = param_configs(6);
num_data_sym = param_configs(7);
% thresholds
corr_threshold = param_configs(10);
pnts_threshold = param_configs(11);

DC = conj(sym_to_data_ang([1],N));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(size(DC_ind,1) == 0)
    return;
end
%% Find list of Potential Preambles from list of Downchirps Detected
pot_pream_ind = [];
c = 1;
for i = 1:size(DC_ind,1)
    if(DC_ind(i,1) - ((num_preamble+num_sync)*N) < 1)
        continue;
    end
    pot_pream_ind(c,:) = DC_ind(i,1) - ((num_preamble + num_sync)*N) : N : DC_ind(i,1)- ((num_sync)*N);
    c = c+1;
end

Upchirp_ind = [];

%% Cross Correlation with a Single UpChirp
temp_wind = [];
for j = 1:size(pot_pream_ind,1)
    if(pot_pream_ind(j,1) - N <= 0)
        continue;
    end
    Data_buffer = [];
    Data_buffer = Data(pot_pream_ind(j,1) - N : pot_pream_ind(j,end)-1 + N);
    temp = [];
    for i = 1:length(Data_buffer) - length(DC)
        temp(i+1) = sum(Data_buffer(i + 1 : i + N).*DC(1:N))...
        / sqrt(sum( Data_buffer(i + 1: i + N) .* conj(Data_buffer(i + 1 : i + N)) ) * ...
        sum( DC(1:N) .* conj(DC(1:N))));
    end
    temp_wind(j,:) = temp;
end

array_stack = {};

% iterate over each Downchirp Detected
for m = 1:size(temp_wind,1) 
    
    n_samp_array = [];
    peak_ind_prev = [];
    for i = 0:floor(length(temp_wind)/N)-1
        % windowing Cross-Correlation Arrays correponsing to each pkt (window length N samples)
        wind = abs(temp_wind(m,i*N + 1 : (i+1) * N));
        peak_ind_curr = get_max(wind,corr_threshold,pnts_threshold);
        if(length(peak_ind_prev) ~= 0 && length(peak_ind_curr) ~= 0)
            for j = 1:length(peak_ind_curr)
                for k = 1:length(peak_ind_prev)
                    % check if combination of any two peaks in consecutive window are N samples apart
                    if(abs(peak_ind_curr(j) == peak_ind_prev(k)))
                        n_samp_array = [n_samp_array  peak_ind_prev(k)+((i-1)*N)+(pot_pream_ind(m,1)-N-1)];
                    end
                % This extracts a list of all peaks that are N samples
                % apart
                end
            end
        end
        peak_ind_prev = peak_ind_curr;
    end
    array_stack{m} = n_samp_array;

end

for m = 1:length(array_stack)
    n_samp_array = [];
    n_samp_array = cell2mat(array_stack(m));
    
    for i = 1:length(n_samp_array)
        c = 0;
        ind_arr = n_samp_array(i) + N : N : n_samp_array(i) + N + ((num_preamble-2)*N);

        for j = 1:length(ind_arr)
            c = c + sum( n_samp_array == ind_arr(j) );
        end
        % Find from the list all the peaks that appear consecutively for
        % more than 6 windows (Upchirp should give 8 peaks, N sampled apart)
        if( c >= 6 )
            if(length(Upchirp_ind) ~= 0)
                if(sum(n_samp_array(i) == Upchirp_ind(:,1)) ~= 1)
                    Upchirp_ind = [Upchirp_ind; [n_samp_array(i) ind_arr]];
                else
                    
                end
            else
                Upchirp_ind = [Upchirp_ind; [n_samp_array(i) ind_arr]];
            end
        end
    end
    
end

% filter Upchirps that are with in 5 samples (same pkt detected multiple times due to peak energy spread)
temp = [];
indices = [zeros(1,num_preamble); Upchirp_ind];
for i = 2:size(indices,1)
    if(length(temp) == 0)
        temp = [temp; indices(i,:)];
    else
        if( min(abs(indices(i) - temp(:,1))) > 5 )
            temp = [temp; indices(i,:)];
        end
    end
end
Upchirp_ind = temp;

end