function [symbols] = CIC_Demod(Pream_ind,Rx_Buffer,Pream_ind_stack,Peak_amp,m)
% CIC Demodulation
% This Function performs the core CIC Demodulation, IT first
% finds the interfering symbol boundaries with in our current demodulation
% window, Finds the potential list of symbols using std dechirping and CIC
% dechirping, filters the list of symbols based on power and frequency
% offset filtering, finds intersection of potential symbols and then makes
% the final decision (refer to the Block Diagram)

% chirp variables
SF = param_configs(1);
N = 2^SF;

% LORA pkt variables
num_preamble = param_configs(4);
num_sync = param_configs(5);
num_DC = param_configs(6);
num_data_sym = param_configs(7);

DC = conj(sym_to_data_ang([1],N));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pream_ind_stack(:,num_preamble + 1) = Pream_ind_stack(:,num_preamble) + N;

% for each Preamble in the Pream_ind_stack, compute exact start and end
% indices for all the data symbols
for i = 1: size(Pream_ind_stack,1)
    frm_st = Pream_ind_stack(i,1) + (num_preamble*N) + (num_DC*N) + (num_sync*N);
    frm_ind(i,:,:) = [((frm_st: N :frm_st+((num_data_sym-1)*N)))' ((frm_st+N-1 : N :frm_st+((num_data_sym)*N)))'];
end

% for the pkt to be demodulated, find indexes for each data symbol
Data_frame_start = Pream_ind(1) + (num_preamble*N) + (num_DC*N) + (num_sync*N);
frame_indices = [((Data_frame_start: N :Data_frame_start+((num_data_sym-1)*N)))' ((Data_frame_start+N-1 : N :Data_frame_start+((num_data_sym)*N)))'];
    
for  k = 1:num_data_sym
    %% Find interfering Symbol Boundaries
    % for current demodulation window, find the chunks of interfering
    % symbols due to collisions by determining index overlaps
    ind = [];
    sym_bnd = [];
    for i = 1:size(frm_ind,1)
        if(i ~= m)
            st = reshape(frm_ind(i,:,1),[],1);
            ed = reshape(frm_ind(i,:,2),[],1);
            sym_bnd = [sym_bnd st(intersect(find(st > frame_indices(k,1)) , find(st < frame_indices(k,2))))];
            sym_bnd = [sym_bnd ed(intersect(find(ed > frame_indices(k,1)) , find(ed < frame_indices(k,2))))];
        end
    end
        
    %% CIC Filtering
    if(frame_indices(k,2) > length(Rx_Buffer))
        symbols(k) = -3;
        continue;
    end
    % standard LoRa dechirping
    data_wind = Rx_Buffer(frame_indices(k,1):frame_indices(k,2)) .* DC;
    data_fft = abs(fft(data_wind));

    % scale the dmeodulation window with appropriate Gaussian to
    % suppress the interfering symbols at window ends
    sigma = 1;
    amp_scale = exp(-(1/(2*(sigma^2)))* linspace(-1,1,N).^2);
    amp_scale = amp_scale./(sqrt(2*pi)*sigma);
    temp_wind = data_wind .* amp_scale;

    sym_bnd = mod(sym_bnd - frame_indices(k,1),N);
    intf_wind = [];
    % n_pnt is the fft Factor - fft(Signal, n_pnt * length(Signal))
    n_pnt = 4;
    for i = 1:length(sym_bnd)
        buff = zeros(2,n_pnt*N);
        buff(1,1:sym_bnd(i) - 1) = temp_wind(1:sym_bnd(i) - 1);
        buff(1,:) = abs(fft(buff(1,:),n_pnt*N))./sqrt(sum(abs(buff(1,:)).^2));
        buff(2,sym_bnd(i):N) = temp_wind(sym_bnd(i):N);
        buff(2,:) = abs(fft(buff(2,:),n_pnt*N))./sqrt(sum(abs(buff(2,:)).^2));
        intf_wind = [intf_wind; buff];
    end
    % CIC's min operation to suppress interfering symbol's Peaks and
    % then finding candidate symbols
    intf_wind_min_fft = min(intf_wind,[],1);
    pot_sym_cic = get_max(intf_wind_min_fft,4*sum(intf_wind_min_fft)/(n_pnt*N),n_pnt*N); % try 3*
    pot_sym_cic = ceil(pot_sym_cic/n_pnt);
    %% Power-Filtering
    % Out of all peaks, find ones with in range of +- 0.5*Preamble Peak
    PwrFctr = 0.5;
    PwrFlr = 4; % may try 3
    up_thresh = (Peak_amp(1) + PwrFctr*Peak_amp(1));
    low_thresh = (Peak_amp(1) - PwrFctr*Peak_amp(1));
    if(low_thresh < (PwrFlr*sum(data_fft)/N)) %1
        low_thresh = (PwrFlr*sum(data_fft)/N);
    end
    pot_sym_pf = get_bounded_max(data_fft,up_thresh,low_thresh);
    %% Filtering Preamble of interfering Packets
    % Filter out Peaks with in current window that are appearinf repeatedly
    % in 3 consecutive windows
    if(~(frame_indices(k,2) + N > length(Rx_Buffer) || frame_indices(k,2) + 2*N > length(Rx_Buffer)))
        data_wind_next_1 = Rx_Buffer(frame_indices(k,1) + N:frame_indices(k,2) + N) .* DC;
        data_wind_prev_1 = Rx_Buffer(frame_indices(k,1) - N:frame_indices(k,2) - N) .* DC;
        data_wind_next_2 = Rx_Buffer(frame_indices(k,1) + 2*N:frame_indices(k,2) + 2*N) .* DC;
        data_wind_prev_2 = Rx_Buffer(frame_indices(k,1) - 2*N:frame_indices(k,2) - 2*N) .* DC;
        temp_next_1 = abs(fft(data_wind_next_1,N));
        temp_prev_1 = abs(fft(data_wind_prev_1,N));
        temp_next_2 = abs(fft(data_wind_next_2,N));
        temp_prev_2 = abs(fft(data_wind_prev_2,N));        
        next_wind_sym_1 = get_max(temp_next_1,4*sum(temp_next_1)/N,N);
        next_wind_sym_2 = get_max(temp_next_2,4*sum(temp_next_2)/N,N);
        prev_wind_sym_1 = get_max(temp_prev_1,4*sum(temp_prev_1)/N,N);
        prev_wind_sym_2 = get_max(temp_prev_2,4*sum(temp_prev_2)/N,N);
        temp = [];
        for i = 1:length(pot_sym_pf)
            if( (sum(pot_sym_pf(i) == prev_wind_sym_1) && sum(pot_sym_pf(i) == next_wind_sym_1))...
                    || (sum(pot_sym_pf(i) == prev_wind_sym_2) && sum(pot_sym_pf(i) == prev_wind_sym_1))...
                    || (sum(pot_sym_pf(i) == next_wind_sym_1) && sum(pot_sym_pf(i) == next_wind_sym_2)) )
            else
                temp = [temp pot_sym_pf(i)];
            end
        end
        pot_sym_pf = temp;
    end
    %%  Freq. Offset Filtering
    % since we have removed Frequency Offset from pkt under consideration 'Pream_ind'
    % and chosen the right downsampled buffer, the true symbol peak should
    % be the most crisp one (either use this or Choir Module(next), results should be almost similar)
    % and the interfering symbols peak may or may not be crisp
    temp = [];
    for i = 1:length(pot_sym_pf)
        if(sum(pot_sym_pf(i) + 1 == pot_sym_pf) || sum(pot_sym_pf(i) - 1 == pot_sym_pf))
        else
            temp = [temp pot_sym_pf(i)];
        end
    end
    pot_sym = temp;
    %%  Choir Module
    % Since the Freq. Offset for the pkt of interest has been compensated,
    % now the Freq. offset of the correct Peak should be closest to 0
%         npnt = 16;
%         data_fft_npnt = abs(fft(data_wind,npnt*N));
%         FO_thresh = 0.25;
%         sym_FO = [];
%         temp = [];
%         for i = 1:length(pot_sym_pf)
%             ind = [];
%             if(pot_sym_pf(i) == 1)
%                 ind = [(N*npnt) - (npnt/2) + 1 : (N*npnt) (((pot_sym_pf(i)-1) * npnt) + 1) + (npnt/2) : N*npnt];
%             else
%                 ind = (((pot_sym_pf(i)-1) * npnt) + 1) - (npnt/2) : (((pot_sym_pf(i)-1) * npnt) + 1) + (npnt/2);
%             end
%             [~,a] = max(data_fft_npnt(ind));
%             sym_FO = [sym_FO abs(a - ((npnt/2)+1))/npnt];
%             if(sym_FO(end) < FO_thresh)
%                 temp = [temp pot_sym_pf(i)];
%             end
%         end
% %         sym(k)
% %         pot_sym_pf - 2
% %         sym_FO
%         pot_sym = temp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         pot_sym = pot_sym_pf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Make the Final Decision
    b = [];
    if(length(sym_bnd) == 0)
        % if there is no symbol colliding with current demod window
        if(length(pot_sym) == 0)
            
            [~,symbols(k)] = max(data_fft);
            
        else

            % choose peak closest in height to Preamble Peak
            dist = abs(data_fft(pot_sym) - (up_thresh + low_thresh)/2);
            [~,b] = min(dist);
            symbols(k) = pot_sym(b);

        end
    else
        % if symbols are colliding with current demod window
        fin_sym = intersect(pot_sym_cic,pot_sym);
        if(length(fin_sym) == 0)
            if(length(pot_sym_cic) == 0 && length(pot_sym) ~= 0)
                
                % choose peak closest in height to Preamble Peak
                dist = abs(data_fft(pot_sym) - (up_thresh + low_thresh)/2);
                [~,b] = min(dist);
                symbols(k) = pot_sym(b);
                
            elseif(length(pot_sym_cic) ~= 0 && length(pot_sym) == 0)
                
                % make decision based on CIC's windows, correct symbol should
                % have lowest std as it appears in all windows
                sdev = std(intf_wind(:,n_pnt.*pot_sym_cic),1);
                [~,b] = min(sdev);
                symbols(k) = pot_sym_cic(b);
                
            elseif(length(pot_sym) == 0 && length(pot_sym_cic) == 0)
                
                [~,symbols(k)] = max(data_fft);
            else
                
                % choose peak closest in height to Preamble Peak
                dist = abs(data_fft(pot_sym) - (up_thresh + low_thresh)/2);
                [~,b] = min(dist);
                symbols(k) = pot_sym(b);

            end
        else
            % if intersection yields some candidates then decide based on partial STFT as following
            
            %%  Stft
            % find stft 2D matrix of follwoing dimensions
            % N frequency (rows)  x  [1 : avg_pnts      N - avg_pnts : N] (columns)
            % i.e. finding Spectrum of first 10 and last 10 time samples (Dont need to compute whole Spectrum)
            avg_pnts = 10;  % # of start and end time samples to average over
            G_wind1 = data_wind;
            Spec = [];
            for i = 0:avg_pnts
                Spec(1:N/2 + i,i+1) = G_wind1(1:N/2 + i);
                Spec(N/2 - (avg_pnts - i):N,i+1 + (avg_pnts + 1)) = G_wind1(N/2 - (avg_pnts - i):N);
            end
            Spec = fft(Spec);
            % the amplitude difference at the start and end of the spectrum
            % for correct symbol should be minimum, (non-interfering symbol's
            % frequency track appears continuosly for all the columns)
            freq_amp = min(abs(Spec(fin_sym,1:(avg_pnts+1))),[],2);
            freq_amp_end = min(abs(Spec(fin_sym,end-avg_pnts:end)),[],2);
            dif = abs(freq_amp - freq_amp_end);
            [~,b] = min(dif);
            symbols(k) = fin_sym(b);

        end
    end
end


end
