close all;
clear all;
clc

%% Loading variables
% chirp variables
SF = param_configs(1);
BW = param_configs(2);
Fs = param_configs(3);
N = 2^SF;
upsampling_factor = Fs / BW;
Ts = 1 / Fs;

% LORA pkt variables
num_preamble = param_configs(4);
num_sync = param_configs(5);
num_DC = param_configs(6);
num_data_sym = param_configs(7);
preamble_sym = 1;
pkt_len = num_preamble + num_sync + num_DC + num_data_sym;
num_samples = pkt_len * N;

% load true symbols
load('sym.mat')                                     

% Generating a Downchirp
DC = conj(sym_to_data_ang([1],N));

%%  Loading File
path = param_configs(14);
fil_nm = param_configs(15);
fi_1 = fopen([path fil_nm]);
x_inter_1 = fread(fi_1,'float32');
fclose(fi_1);

% parse complex data
x_1 = x_inter_1(1:2:end) + 1i * x_inter_1(2:2:end);	% Create complex values
x_1 = x_1.';

% file Duration
t = [0:length(x_1)-1] / Fs;

x_1 = x_1(1:floor(length(x_1) / upsampling_factor) * upsampling_factor);
x_1_dnsamp = x_1(1:upsampling_factor:end);
file_dur = (length(x_1)/Fs);

% plotting downsampled signal
plot(linspace(0,file_dur,length(x_1_dnsamp)),real(x_1_dnsamp),'linewidth',0.75)
xlabel('time/sec')
ylabel('real magnitude')

%%  Active Sessions Detection using Dechirping Windows
uplink_wind = active_sess_dechirp(x_1);		% uplink_wind contains the [start,  end] indices of each active session detected
disp(['Detected ' num2str(size(uplink_wind,1)) ' active sessions'])
%%
demod_sym_stack = [];
Peaks = [];
%%
for m = 1:size(uplink_wind,1)
    disp(' ')
    disp(['Active Session no. ' num2str(m)])

%%      DC correlations to find LoRa pkts out of collision
    temp_buff = [];
    temp_buff = x_1(uplink_wind(m,1) : uplink_wind(m,2));
    temp_buff = temp_buff(:,1:floor(size(temp_buff,2) / upsampling_factor) * upsampling_factor);

    DC_ind = DC_location_correlation(temp_buff(1:upsampling_factor:end));
    disp(['Found ' num2str(size(DC_ind,1)) ' Downchirps in current Active session'])
    
    if(size(DC_ind,1) == 0)
        continue;
    end

%%      UC correlation to filter false positives and frequency offset & Packets' SNR estimation
    % All possible downsampled Buffers with different starting sample for downsampling
    Data_freq_off = [];
    Rx_Buff_dnsamp = [];
    for i = 1:upsampling_factor
        Rx_Buff_dnsamp(i,:) = temp_buff(i:upsampling_factor:end);
    end
    [Upchirp_ind] = UC_location_corr_DC_based(temp_buff(1:upsampling_factor:end),DC_ind);
    
    if(size(Upchirp_ind,1) == 0)
        continue;
    end
    
    % for each preamble detected, choose the correct downsampled buffer with any frequency offset
    % compensated and determine preamble peak heights to be used later for power-filtering
    [Data_freq_off, Peak, Upchirp_ind,FFO] = dnsamp_buff(Rx_Buff_dnsamp,Upchirp_ind);

    if(size(Upchirp_ind,1) == 0)
        continue;
    end
    
    % Filter False Positives based on 2-SYNC Words detected 
    [Preamble_ind, bin_offsets, Data_out, Peak_amp,FFO] = filter_false_postives(Data_freq_off,Upchirp_ind,Peak,FFO);
%%  filter preambles that are with in 5 samples (same pkt detected twice due to Correlation peak energy spread)
    temp = [];
    temp_data = [];
    temp_peaks = [];
    indices = [zeros(1,num_preamble); Preamble_ind];
    Data = [zeros(1,size(Data_out,2)); Data_out];
    peaks = [zeros(1,size(Peak_amp,2)); Peak_amp];
    clear Peak_amp
    for i = 2:size(indices,1)
        if(length(temp) == 0)
            temp = [temp; indices(i,:)];
            temp_data = [temp_data; Data(i,:)];
            temp_peaks = [temp_peaks; peaks(i,:)];
        else
            if( min(abs(indices(i) - temp(:,1))) > 5 )
                temp = [temp; indices(i,:)];
                temp_data = [temp_data; Data(i,:)];
                temp_peaks = [temp_peaks; peaks(i,:)];
            end
        end
    end
    Pream_ind = temp;
    Data_out = temp_data;
    Peak_amp = temp_peaks;
    
    disp(['Found ' num2str(size(Pream_ind,1)) ' Preambles in current Active session'])

    %%  Data Demodulation using CIC
    demod_sym = [];
    for j =1:size(Pream_ind,1)
        disp(['demodulating ' num2str(j) 'th pkt in current active session'])
        [demod_sym(j,:)] = CIC_Demod(Pream_ind(j,:),Data_out(j,:),Pream_ind,Peak_amp(j,:),j);
        demod_sym(j,:) = mod(demod_sym(j,:) + bin_offsets(j) - 2,N);
    end
    demod_sym_stack = [demod_sym_stack; demod_sym];
    Peaks = [Peaks; Peak_amp];
end

%%  Calculating SER for the File
ser = sum(sum(repmat(sym,size(demod_sym_stack(4:end,:),1),1) ~= demod_sym_stack(4:end,:))) / prod(size(demod_sym_stack(4:end,:)));
disp('******************************************')
disp(['Symbol Error Rate for this File = ' num2str(ser)])
disp('*****************FINISHED*****************')
%%  Calculating Symbol Throughput for the File
ser = sum(sum(repmat(sym,size(demod_sym_stack(4:end,:),1),1) == demod_sym_stack(4:end,:))) / 60;    % 60 = duration of each aggregate Rate session
disp('******************************************')
disp(['Symbol Throughput for this file = ' num2str(ser)])
disp('*****************FINISHED*****************')