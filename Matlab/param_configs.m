function par = param_configs(id)

    % LoRa PHY transmitting parameters
    LORA_SF = 8;		% LoRa spreading factor - 7, 8, 9, 10, 11, 12
    LORA_BW = 250e3;		% LoRa bandwidth - 125e3, 250e3, 500e3
    
    % Receiving device parameters
    Fs = 2e6;       		% recerver's sampling rate 
    num_preamble = 8;		% num of Preamble Base Upchirps in a Lora Pkt

    num_sync = 2;
    num_DC = 2.25;
    num_data_sym = 28;
    
    DC_corr_threshold = 0.2;        % Value to be calibrated based on Correlation plot's noise floor
    DC_corr_pnts_threshold = 40;
    
    UC_corr_threshold = 0.1;        % Value to be calibrated based on Correlation plot's noise floor
    UC_corr_pnts_threshold = 40;
    SYNC1 = 8;
    SYNC2 = 16;
    
    path = 'C:\Matlab_lt\CIC\low snr';            % Add path to the file
    fil_nm = '\60sec_perpkt_lowsnr';                    % File name
    
    
    switch(id)
        case 1,
            par = LORA_SF;
        case 2,
            par = LORA_BW;
        case 3,
            par = Fs;
        case 4,
            par = num_preamble;
        case 5,
            par = num_sync;
        case 6,
            par = num_DC;
        case 7,
            par = num_data_sym;
        case 8,
            par = DC_corr_threshold;
        case 9,
            par = DC_corr_pnts_threshold;
        case 10,
            par = UC_corr_threshold;
        case 11,
            par = UC_corr_pnts_threshold;
        case 12,
            par = SYNC1;
        case 13,
            par = SYNC2;
        case 14,
            par = path;
        case 15,
            par = fil_nm;

        otherwise,
    end
end
