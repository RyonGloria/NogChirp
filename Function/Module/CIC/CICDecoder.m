classdef CICDecoder
    properties
        loraSet;
        binRecord;
        detectCount;
    end
    
    methods
        function obj = CICDecoder(loraSet)
            obj.loraSet = loraSet;
        end
        
        function obj = decode(obj, x_1)
            %% Loading variables
            % chirp variables
            obj.binRecord = cell(1, 0);
            SF = obj.param_configs(1);
            BW = obj.param_configs(2);
            Fs = obj.param_configs(3);
            N = 2^SF;
            upsampling_factor = Fs/BW;
            Ts = 1/Fs;
            
            % LORA pkt variables
            num_preamble = obj.param_configs(4);
            num_sync = obj.param_configs(5);
            num_DC = obj.param_configs(6);
            num_data_sym = obj.param_configs(7);
            preamble_sym = 1;
            pkt_len = num_preamble + num_sync + num_DC + num_data_sym;
            num_samples = pkt_len * N;
            
            % Generating a Downchirp
            DC = conj(obj.sym_to_data_ang([1],N));
            
            %%  Loading File
            % path = param_configs(14);
            % fil_nm = param_configs(15);
            % path = 'F:\Pyramid_temp\SF9BW125\RTL4';            % Add path to the file
            % fil_nm = '\T17_26_30_SF9_BW125000.sigmf-data';
            % fi_1 = fopen([path fil_nm]);
            % x_inter_1 = fread(fi_1,'float32');
            % fclose(fi_1);
            
            % parse complex data
            % x_1 = x_inter_1(1:2:end) + 1i*x_inter_1(2:2:end);   % Read Complex values
            % x_1 = x_1.';
            % x_1 = x_1 + 0.5*circshift(x_1, N*upsampling_factor*20.75);
            % signalLength =
            % x_1 = [x_1, zeros(1, N*upsampling_factor*500)];
            % x_1 = [x_1, x_1];
            % x_1 = x_1(1:N*upsampling_factor*40);
            
            % file Duration
            t = [0:length(x_1)-1]/Fs;
            
            x_1 = x_1(1:floor(length(x_1)/upsampling_factor)*upsampling_factor);
            x_1_dnsamp = x_1(1:upsampling_factor:end);
            file_dur = (length(x_1)/Fs);
            
            %%  Active Sessions Detection using Dechirping Windows
            uplink_wind = obj.active_sess_dechirp(x_1);             % uplink_wind contains the [start,  end] indices of each active session detected
            uplink_wind = obj.active_sess_split(uplink_wind, 10 * num_samples * upsampling_factor, 2.5 * upsampling_factor * N);    % split up sessions longer than 10 times packet length
            % disp(['Detected ' num2str(size(uplink_wind,1)) ' active sessions'])
            %%
            demod_sym_stack = [];
            Peaks = [];
            %%
            for m = 1:size(uplink_wind,1)
                % disp(' ')
                % disp(['Active Session no. ' num2str(m)])
                
                %%      DC correlations to find LoRa pkts out of collision
                
                temp_buff = [];
                temp_buff = x_1(uplink_wind(m,1) : uplink_wind(m,2));
                temp_buff = temp_buff(:,1:floor(size(temp_buff,2)/upsampling_factor)*upsampling_factor);
                
                DC_ind = obj.DC_location_correlation(temp_buff(1:upsampling_factor:end));
                % disp(['Found ' num2str(size(DC_ind,1)) ' Downchirps in current Active session'])
                
                if(size(DC_ind,1) == 0)
                    continue;
                end
                if((DC_ind(end,1)*upsampling_factor + ((num_DC + num_data_sym)*N*upsampling_factor)) > length(temp_buff))
                    %   if a data portion of a packet is split then account for the
                    %   length difference
                    ex_samp = (DC_ind(end,1)*upsampling_factor + ((num_DC + num_data_sym)*N*upsampling_factor)) - length(temp_buff);
                    temp_buff = x_1(uplink_wind(m,1) : uplink_wind(m,2) + ex_samp);
                    temp_buff = temp_buff(:,1:floor(size(temp_buff,2)/upsampling_factor)*upsampling_factor);
                end
                
                
                %%      UC correlation to filter false positives and frequency offset & Packets' SNR estimation
                % All possible downsampled Buffers with different starting sample for downsampling
                Data_freq_off = [];
                Rx_Buff_dnsamp = [];
                for i = 1:upsampling_factor
                    Rx_Buff_dnsamp(i,:) = temp_buff(i:upsampling_factor:end);
                end
                [Upchirp_ind] = obj.UC_location_corr_DC_based(temp_buff(1:upsampling_factor:end),DC_ind);
                
                if(size(Upchirp_ind,1) == 0)
                    continue;
                end
                
                % for each Preamble detected, Choosing the correct downsampled buffer with any frequency offset
                % been compensated and determining Preamble Peak heights to be used later for Power filtering
                [Data_freq_off, Peak, Upchirp_ind,FFO] = obj.dnsamp_buff(Rx_Buff_dnsamp,Upchirp_ind);
                
                if(size(Upchirp_ind,1) == 0)
                    continue;
                end
                
                % Filter False Positives based on 2-SYNC Words detected
                [Preamble_ind, bin_offsets, Data_out, Peak_amp,FFO] = obj.filter_false_postives(Data_freq_off,Upchirp_ind,Peak,FFO);
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
                
                % disp(['Found ' num2str(size(Pream_ind,1)) ' Preambles in current Active session'])
                %%  Data Demodulation using CIC
                demod_sym = [];
                for j =1:size(Pream_ind,1)
%                     disp(['demodulating ' num2str(j) 'th pkt in current active session'])
                    [demod_sym(j,:)] = obj.CIC_Demod(Pream_ind(j,:),Data_out(j,:),Pream_ind,Peak_amp(j,:),j);
                    demod_sym(j,:) = mod(demod_sym(j,:) + bin_offsets(j) - 2,N);
                end
                demod_sym_stack = [demod_sym_stack; demod_sym];
                Peaks = [Peaks; Peak_amp];
            end

            % 处理bin值
            obj.binRecord = cell(1, size(demod_sym_stack, 1));
            for index = 1:size(demod_sym_stack, 1)
                obj.binRecord{index} = mod(demod_sym_stack(index, :)+2, obj.loraSet.fft_x);
            end
            
        end

        function obj = detectPass(obj, x_1)
            obj.detectCount = 0;
            %% Loading variables
            % chirp variables
            obj.binRecord = cell(1, 0);
            SF = obj.param_configs(1);
            BW = obj.param_configs(2);
            Fs = obj.param_configs(3);
            N = 2^SF;
            upsampling_factor = Fs/BW;
            Ts = 1/Fs;
            
            % LORA pkt variables
            num_preamble = obj.param_configs(4);
            num_sync = obj.param_configs(5);
            num_DC = obj.param_configs(6);
            num_data_sym = obj.param_configs(7);
            preamble_sym = 1;
            pkt_len = num_preamble + num_sync + num_DC + num_data_sym;
            num_samples = pkt_len * N;
            
            % Generating a Downchirp
            DC = conj(obj.sym_to_data_ang([1],N));
            
            % file Duration
            t = [0:length(x_1)-1]/Fs;
            
            x_1 = x_1(1:floor(length(x_1)/upsampling_factor)*upsampling_factor);
            x_1_dnsamp = x_1(1:upsampling_factor:end);
            file_dur = (length(x_1)/Fs);
            
            %%  Active Sessions Detection using Dechirping Windows
            uplink_wind = obj.active_sess_dechirp(x_1);             % uplink_wind contains the [start,  end] indices of each active session detected
            uplink_wind = obj.active_sess_split(uplink_wind, 10 * num_samples * upsampling_factor, 2.5 * upsampling_factor * N);    % split up sessions longer than 10 times packet length
            % disp(['Detected ' num2str(size(uplink_wind,1)) ' active sessions'])
            %%
            demod_sym_stack = [];
            Peaks = [];
            %%
            for m = 1:size(uplink_wind,1)
                % disp(' ')
                % disp(['Active Session no. ' num2str(m)])
                
                %%      DC correlations to find LoRa pkts out of collision
                
                temp_buff = [];
                temp_buff = x_1(uplink_wind(m,1) : uplink_wind(m,2));
                temp_buff = temp_buff(:,1:floor(size(temp_buff,2)/upsampling_factor)*upsampling_factor);
                
                DC_ind = obj.DC_location_correlation(temp_buff(1:upsampling_factor:end));
                % disp(['Found ' num2str(size(DC_ind,1)) ' Downchirps in current Active session'])
                
                if(size(DC_ind,1) == 0)
                    continue;
                end
                if((DC_ind(end,1)*upsampling_factor + ((num_DC + num_data_sym)*N*upsampling_factor)) > length(temp_buff))
                    %   if a data portion of a packet is split then account for the
                    %   length difference
                    ex_samp = (DC_ind(end,1)*upsampling_factor + ((num_DC + num_data_sym)*N*upsampling_factor)) - length(temp_buff);
                    temp_buff = x_1(uplink_wind(m,1) : uplink_wind(m,2) + ex_samp);
                    temp_buff = temp_buff(:,1:floor(size(temp_buff,2)/upsampling_factor)*upsampling_factor);
                end
                
                
                %%      UC correlation to filter false positives and frequency offset & Packets' SNR estimation
                % All possible downsampled Buffers with different starting sample for downsampling
                Data_freq_off = [];
                Rx_Buff_dnsamp = [];
                for i = 1:upsampling_factor
                    Rx_Buff_dnsamp(i,:) = temp_buff(i:upsampling_factor:end);
                end
                [Upchirp_ind] = obj.UC_location_corr_DC_based(temp_buff(1:upsampling_factor:end),DC_ind);
                
                if(size(Upchirp_ind,1) == 0)
                    continue;
                end
                
                % for each Preamble detected, Choosing the correct downsampled buffer with any frequency offset
                % been compensated and determining Preamble Peak heights to be used later for Power filtering
                [Data_freq_off, Peak, Upchirp_ind,FFO] = obj.dnsamp_buff(Rx_Buff_dnsamp,Upchirp_ind);
                
                if(size(Upchirp_ind,1) == 0)
                    continue;
                end
                
                % Filter False Positives based on 2-SYNC Words detected
                [Preamble_ind, bin_offsets, Data_out, Peak_amp,FFO] = obj.filter_false_postives(Data_freq_off,Upchirp_ind,Peak,FFO);
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
                obj.detectCount = obj.detectCount + size(Pream_ind, 1);
            end
            
        end
        
        
        function [uplink_wind] = active_sess_dechirp(obj, x_1)
            SF = obj.param_configs(1);
            BW = obj.param_configs(2);
            Fs = obj.param_configs(3);
            N = 2^SF;
            upsampling_factor = Fs/BW;
            
            DC = conj(obj.sym_to_data_ang([1],N));
            DC_fft = fft(DC);
            DC_upsamp =(ifft([DC_fft(1:N/2) zeros(1,(upsampling_factor-1)*N) DC_fft(N/2 + 1:N)]));
            
            peak_gain = [];
            uplink_wind = [];
            n = [];
            p = [];
            last_wind = 0;
            win_jump_factor = 3;
            front_buf = 6*win_jump_factor;
            back_buf = 3*win_jump_factor;
            win_jump = floor(N*upsampling_factor/win_jump_factor);
            mov_thresh_wind = 1000*win_jump_factor;
            mov_thresh = 0;
            mov_thresh_rec  = [];
            for i = 1:floor(length(x_1)/win_jump) - win_jump_factor%length(DC))
                wind = x_1((i-1)*win_jump+1 : (i-1)*win_jump+(N*upsampling_factor));
                wind_fft = abs(fft(wind.*DC_upsamp));
                wind_fft = wind_fft([1:N/2 (N/2 + (upsampling_factor-1)*N)+1:(upsampling_factor)*N]);
                noise_floor = mean(wind_fft);
                n = [n noise_floor];
                fft_peak = max(wind_fft);
                p = [p fft_peak];
                peak_gain = [peak_gain 10*log10(fft_peak/noise_floor)];
                if(i > mov_thresh_wind)
                    mov_thresh = 1.3*mean(peak_gain(end - mov_thresh_wind + 1: end));
                    if(mov_thresh > 6)
                        mov_thresh = 6;
                    end
                else
                    mov_thresh = 1.3*mean(peak_gain);       % 1.3 is adjustable based on noise floor
                    if(mov_thresh > 6)
                        mov_thresh = 6;
                    end
                end
                mov_thresh_rec = [mov_thresh_rec mov_thresh];
                
                if(peak_gain(end) >= mov_thresh)
                    if(i > last_wind)
                        if(i-back_buf < 1)
                            uplink_wind = [uplink_wind; 1 i+front_buf];
                        else
                            uplink_wind = [uplink_wind; i-back_buf i+front_buf];
                        end
                        last_wind = uplink_wind(end,2);
                    elseif(i <= last_wind)
                        uplink_wind(end,2) = i + front_buf;
                        last_wind = uplink_wind(end,2);
                    end
                end
            end
            uplink_wind = uplink_wind(find(uplink_wind(:,2)-uplink_wind(:,1) ~= front_buf + back_buf),:);
            uplink_wind = uplink_wind(find(uplink_wind(:,2)-uplink_wind(:,1) ~= (front_buf + back_buf + 1)),:);
            uplink_wind = uplink_wind(find(uplink_wind(:,2)-uplink_wind(:,1) ~= (front_buf + back_buf - 1)),:);
            temp_link = uplink_wind;
            uplink_wind = uplink_wind.*(win_jump);
            
            if(uplink_wind(end,2) > length(x_1))
                uplink_wind(end,2) = length(x_1);
            end
        end
        
        function [symbols] = CIC_Demod(obj, Pream_ind,Rx_Buffer,Pream_ind_stack,Peak_amp,m)
            % chirp variables
            SF = obj.param_configs(1);
            N = 2^SF;
            
            % LORA pkt variables
            num_preamble = obj.param_configs(4);
            num_sync = obj.param_configs(5);
            num_DC = obj.param_configs(6);
            num_data_sym = obj.param_configs(7);
            
            DC = conj(obj.sym_to_data_ang([1],N));
            
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
                pot_sym_cic = obj.get_max(intf_wind_min_fft,4*sum(intf_wind_min_fft)/(n_pnt*N),n_pnt*N); % try 3*
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
                pot_sym_pf = obj.get_bounded_max(data_fft,up_thresh,low_thresh);
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
                    next_wind_sym_1 = obj.get_max(temp_next_1,4*sum(temp_next_1)/N,N);
                    next_wind_sym_2 = obj.get_max(temp_next_2,4*sum(temp_next_2)/N,N);
                    prev_wind_sym_1 = obj.get_max(temp_prev_1,4*sum(temp_prev_1)/N,N);
                    prev_wind_sym_2 = obj.get_max(temp_prev_2,4*sum(temp_prev_2)/N,N);
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
        
        function [windows] = active_sess_split(obj, windows, max_window_length, window_overlap)
            i = 0;
            while i < size(windows, 1)
                i = i + 1;
                if(abs(windows(i, 2) - windows(i, 1)) > max_window_length)  % need to split window in half
                    len = abs(windows(i, 2) - windows(i, 1)) / 2;               % get bisected length
                    win1_end = ceil(windows(i, 1) + len + (window_overlap / 2));        % add overlap to each half
                    win2_start = floor(windows(i, 2) - len - (window_overlap / 2));
                    
                    windows = [windows(1:i - 1,:); [windows(i,1), win1_end]; windows(i:end,:)]; % bisect the window
                    windows(i+1, 1) = win2_start;
                    i = i - 1;  % see if the first bisected half is within the max window length
                end
            end
        end
        
        function [Downchirp_ind] = DC_location_correlation(obj, Rx_Buffer)
            %DC_LOCATION Summary of this function goes here
            % this function runs the cross-correlation of Rx_Buffer with a single
            % Downchirp and outputs the indices of any Downchirp detected based on
            % 2 correlation peaks that are N samples apart
            
            % loading variables
            SF = obj.param_configs(1);
            BW = obj.param_configs(2);
            Fs = obj.param_configs(3);
            N = 2^SF;
            num_DC = obj.param_configs(6);
            % thresholds
            corr_threshold = obj.param_configs(8);      % Threshold above which we extract all Correlation peaks
            pnts_threshold = obj.param_configs(9);      % Max. # of peaks to extract from Corrrelation Plot
            
            DC = conj(obj.sym_to_data_ang([1],N));
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
            
            n_samp_array = [];
            peak_ind_prev = [];
            for i = 0:floor(length(Cross_Corr)/N)-1
                % windowing Cross-Correlation (window length N samples)
                wind = abs(Cross_Corr(i*N + 1 : (i+1) * N));
                % Extract Multiple Correlation Peaks
                peak_ind_curr = obj.get_max(wind,corr_threshold,pnts_threshold);
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
        
        function [Data_buff peak_amp Up_ind FFO] = dnsamp_buff(obj, Data_stack,Upchirp_ind)
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
            SF = obj.param_configs(1);
            BW = obj.param_configs(2);
            Fs = obj.param_configs(3);
            N = 2^SF;
            num_preamble = obj.param_configs(4);
            num_sync = obj.param_configs(5);
            num_DC = obj.param_configs(6);
            
            DC = conj(obj.sym_to_data_ang([1],N));
            
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
                    Spec = obj.stft_v1(Data_freq_off(m,Upchirp_ind(k,1) - N:Upchirp_ind(k,end) + N - 1 - N),N,DC(1:N),0,0);
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
        
        function [Preamble_ind, bin_offsets, Data_out, Peak_amp, ffo] = filter_false_postives(obj, Data_stack,Upchirp_ind,Peak,FFO)
            %FILTER_FALSE_POSTIVES Summary of this function goes here
            % This function takes in detected Possible Preambles and there relevant
            % data_stacks and then looks at points of concern for the presence of 2
            % SYNC-WORDS. This final filter removes false Positives with a very high
            % certanity
            
            % load parameters
            SF = obj.param_configs(1);
            BW = obj.param_configs(2);
            Fs = obj.param_configs(3);
            N = 2^SF;
            num_preamble = obj.param_configs(4);
            num_sync = obj.param_configs(5);
            num_DC = obj.param_configs(6);
            S1 = obj.param_configs(12);
            S2 = obj.param_configs(13);
            
            DC = conj(obj.sym_to_data_ang(ones(1,num_preamble),N));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ffo = [];
            Preamble_ind = [];
            bin_offsets = [];
            Data_out = [];
            Peak_amp = [];
            % row_ind contains the Frequency Bins around bin 1 where a
            % Preamble Peak can lie, assumption (-6*BW/2^SF <= Freq_off <= 6*BW/2^SF)
            row_ind = [N-4:N+1 1:6];
            for m = 1:size(Upchirp_ind,1)
                % extract 8 Preamble long window
                data_wind = Data_stack(m,Upchirp_ind(m,1) : Upchirp_ind(m,1) + ((num_preamble )*N) -1);
                close all
                % Compute STFT to accurately find the Preamble Frequency Bin and any
                % bin offset (if any left)
                [Spec] = obj.stft_v2(data_wind.*DC,N);
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
                sync_threshold_up = Peak(m,1) + 0.5*Peak(m,1);
                sync_threshold_low = Peak(m,1) - 0.5*Peak(m,1);
                
                sync_word1 = abs(fft(sync_wind(1:N).*DC(1:N)));
                sync_word2 = abs(fft(sync_wind(N+1:end).*DC(1:N)));
                
                if(sync_threshold_low < (2*sum(sync_word1)/N))
                    sync_threshold_low = (2*sum(sync_word1)/N);
                elseif( sync_threshold_low < (2*sum(sync_word2)/N))
                    sync_threshold_low = (2*sum(sync_word2)/N);
                end
                
                % Extract Peaks qualifying Peak Thresholds
                syn1_pnts = obj.get_bounded_max(sync_word1,sync_threshold_up,sync_threshold_low);
                syn2_pnts = obj.get_bounded_max(sync_word2,sync_threshold_up,sync_threshold_low);
                % Check If both SYNC WORDS are present
                if(sum(syn1_pnts == sync1_ind) && sum(syn2_pnts == sync2_ind))
                    Preamble_ind = [Preamble_ind; Upchirp_ind(m,:)];
                    if(pream_peak_ind(m) < N/2)
                        bin_offsets = [bin_offsets 1 + (-mod(pream_peak_ind(m),N))];
                    else
                        bin_offsets = [bin_offsets mod(N+2 - pream_peak_ind(m),N)];
                    end
                    Data_out = [Data_out; Data_stack(m,:)];
                    Peak_amp = [Peak_amp; Peak(m,:)];
                    ffo = [ffo; FFO(m)];
                end
            end
            
            
        end
        
        function [pnts] = get_bounded_max(obj, arr,up_thresh,low_thresh)
            %GET_BOUNDED_MAX, this function returns all points in arr that are >
            %low_thresh and < up_thresh
            
            pnts_up = find(arr < up_thresh);
            pnts_low = find(arr > low_thresh);
            pnts = intersect(pnts_up,pnts_low);
        end
        
        function [out] = get_max(obj, arr,threshold,num_pnts)
            %GET_MAX, this function returns uptil num_pnts many maximums above
            % 'threshold'
            
            out = [];
            for i = 1:num_pnts
                [a b] = max(arr);
                if(a < threshold)
                    return;
                else
                    out = [out b];
                    arr(b) = 0;
                end
            end
            
        end
        
        function par = param_configs(obj, id)
            
            % LoRa PHY transmitting parameters
            LORA_SF = obj.loraSet.sf;%12;            % LoRa spreading factor
            LORA_BW = obj.loraSet.bw;%125e3;        % LoRa bandwidth
            
            % Receiving device parameters
            Fs = obj.loraSet.sample_rate;%125e3*8;  % recerver's sampling rate
            num_preamble = obj.loraSet.Preamble_length;       % num of Preamble Base Upchirps in a Lora Pkt
            num_sync = 2;
            num_DC = 2.25;
            num_data_sym = obj.loraSet.payloadNum;
            
            DC_corr_threshold = 0.2;        % Value to be calibrated based on Correlation plot's noise floor
            DC_corr_pnts_threshold = 40;
            
            UC_corr_threshold = 0.1;        % Value to be calibrated based on Correlation plot's noise floor
            UC_corr_pnts_threshold = 40;
            SYNC1 = 8;
            SYNC2 = 16;
            
            %     path = 'F:\Pyramid_temp\SF9BW125\RTL4';            % Add path to the file
            %     fil_nm = '\T17_26_30_SF9_BW125000.sigmf-data';                    % File name
            
            
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
        
        function [Spec] = stft_v1(obj, Rx_Buffer,N,DC,upsamp,dis)
            %STFT
            % This function produces a spectrogram of LoRa signal using Dechirping
            % operation to get the best frequency Resolution
            Spec = zeros(N,length(Rx_Buffer));
            buff = [Rx_Buffer zeros(1,N-1)];
            if(~upsamp)
                for i = 1:length(Rx_Buffer)
                    Spec(:,i) = circshift(abs(fft(buff(i:i+N-1).*DC))./sqrt(N),-(i-1));%\
                end
                if(dis == 1)
                    %             spec_plot(abs(Spec),N,0,0,0)
                end
            end
        end
        
        function [Spec] = stft_v2(obj, Buffer,f_n)
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
        
        function [data] = sym_to_data_ang(obj, symbol,N)
            %SYM_TO_DATA returns an N-sample upchirp of data symbol
            
            data = [];
            accumulator = 0;
            
            for j = symbol
                phase = -pi + ((j-1)*(2*pi/(N)));
                temp = [];
                for i = 1:N
                    accumulator = accumulator + phase;
                    polar_radius = 1;
                    
                    [x, y] = pol2cart(accumulator, polar_radius);
                    
                    temp(i) = complex(x, y);
                    
                    phase = phase + (2*pi/(N));
                end
                data = [data temp];
            end
        end
        
        function symbol_reformat(obj, demod_sym_stack)
            %% write symbols to txt file, to be fed to RPP0
            % Each Symbol appears on a column and symbols of different packets are separated
            % with -1 identifier
            fileID_1 = fopen('symbols.txt','w');
            demod_sym_stack(find(demod_sym_stack < 0)) = 0;
            for i = 1:size(demod_sym_stack,1)
                fprintf(fileID_1,'%d\n',-1);
                fprintf(fileID_1,'%d\n',demod_sym_stack(i,:));
            end
            fclose(fileID_1);
            
        end
        
        function [Upchirp_ind] = UC_location_corr_DC_based(obj, Data,DC_ind)
            %UC_location_corr_DC_based
            % this function takes in the Data buffer, the Downchirp indices from
            % previous function and looks locally for presence of an Upchirp in order
            % to increase the certainity of the presence of a pkt
            
            % chirp variables
            SF = obj.param_configs(1);
            BW = obj.param_configs(2);
            Fs = obj.param_configs(3);
            N = 2^SF;
            upsampling_factor = Fs/BW;
            
            % LORA pkt variables
            num_preamble = obj.param_configs(4);
            num_sync = obj.param_configs(5);
            num_DC = obj.param_configs(6);
            num_data_sym = obj.param_configs(7);
            % thresholds
            corr_threshold = obj.param_configs(10);
            pnts_threshold = obj.param_configs(11);
            
            DC = conj(obj.sym_to_data_ang([1],N));
            
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
                    peak_ind_curr = obj.get_max(wind,corr_threshold,pnts_threshold);
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
        
        
    end
    
end