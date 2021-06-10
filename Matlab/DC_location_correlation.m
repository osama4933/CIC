function [Downchirp_ind] = DC_location_correlation(Rx_Buffer)
%DC_LOCATION Summary of this function goes here
% this function runs the cross-correlation of Rx_Buffer with a single
% Downchirp and outputs the indices of any Downchirp detected based on
% 2 correlation peaks that are N samples apart

% loading variables
SF = param_configs(1);
BW = param_configs(2);
Fs = param_configs(3);
N = 2^SF;
num_DC = param_configs(6);
% thresholds
corr_threshold = param_configs(8);      % Threshold above which we extract all Correlation peaks
pnts_threshold = param_configs(9);      % Max. # of peaks to extract from Corrrelation Plot

DC = conj(sym_to_data_ang([1],N));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cross Correlation with a Single downchirp
Downchirp_ind = [];
for i = 1:length(Rx_Buffer) - length(DC) - 1
    Cross_Corr(i) = sum(Rx_Buffer( i : i + (N) - 1) .* conj(DC))...
            / sqrt(sum( Rx_Buffer( i : i + (N) - 1) .* conj(Rx_Buffer( i : i + (N) - 1)) ) * ...
            sum( DC .* conj(DC)));
end
Cross_Corr = Cross_Corr(isfinite(Cross_Corr));
corr_threshold =  4*sum(abs(Cross_Corr))/length(Cross_Corr);

%%  Optional Cross-Correlation Plot
% figure
% plot(abs(Cross_Corr))
% % [~,Downchirp_ind] = max(Cross_Corr);
% hold on
% plot(corr_threshold.*ones(1,length(Cross_Corr)))
% set(gca,'linewidth',1.5,'fontsize',25,'fontname','Times New Roman');
% title('Correlation with single Downchirp','FontSize',30);
% xlabel('Samples','FontSize',30);
% ylabel('Amp.','FontSize',30);
% ylim([0 1])

% keyboard

n_samp_array = [];
peak_ind_prev = [];
for i = 0:floor(length(Cross_Corr)/N)-1
    % windowing Cross-Correlation (window length N samples)
    wind = abs(Cross_Corr(i*N + 1 : (i+1) * N));                            
    % Extract Multiple Correlation Peaks
    peak_ind_curr = get_max(wind,corr_threshold,pnts_threshold);          
    if(length(peak_ind_prev) ~= 0 && length(peak_ind_curr) ~= 0)
        for j = 1:length(peak_ind_curr)
            for k = 1:length(peak_ind_prev)
                % check if combination of any two peaks in consecutive window are N samples apart
                if(peak_ind_curr(j) == peak_ind_prev(k))                    
                    n_samp_array = [n_samp_array  peak_ind_prev(k)+((i-1)*N) peak_ind_curr(j)+((i)*N)];
                end
                % This extracts a list of all peaks that are N samples
                % apart
            end
        end 
    end
    peak_ind_prev = peak_ind_curr;
end

for i = 1:length(n_samp_array)
    c = 0;
    ind_arr = n_samp_array(i) : N : n_samp_array(i) + (N);
    
    for j = 1:length(ind_arr)
        c = c + sum( n_samp_array == ind_arr(j) );
    end
    % Find from the list all the peaks that appear consecutively for
    % more than 2 windows (Downchirp should give 3 peaks, N sampled apart)
    if( c >= 2 )
        Downchirp_ind = [Downchirp_ind; [ind_arr]];
    end
end


% filter Downchirps that are with in 3 samples (same pkt detected twice due to peak energy spread)
temp = [];
indices = [zeros(1,floor(num_DC)); Downchirp_ind];
for i = 2:size(indices,1)
    
    if(isempty(temp))
        temp = [temp; indices(i,:)];
    else
        if( min(abs(indices(i) - temp(:,1))) > 3 )
            temp = [temp; indices(i,:)];
        end
    end
end
Downchirp_ind = temp;


end

