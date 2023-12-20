% 补零计算lora的cfo和time_offset
function [cfo, windows_offset] = get_cfo_winoff_bybin(lora_set, upchirp_peak, downchirp_peak)
    % 计算主峰的CFO(需要补零操作)
    % 对Preamble阶段的FFT峰值进行排序，得到前filter的峰
    zeropadding_size = lora_set.factor;                   % 设置补零的数量，这里的decim表示，补上decim-1倍窗口的零，计算FFT时一共是decim倍的窗口（decim+1）
    d_sf = lora_set.sf;
    d_bw = lora_set.bw;
    fft_x = lora_set.fft_x;
    fft_x_zeropadding = fft_x*zeropadding_size*2^(10-d_sf);
    
    % 计算CFO和窗口偏移量
    if upchirp_peak + downchirp_peak < fft_x_zeropadding*0.5
        cfo_bin = upchirp_peak + downchirp_peak - 2;
        cfo = -cfo_bin/2/fft_x_zeropadding * d_bw;
        windows_offset = (downchirp_peak - upchirp_peak) / 2^(11-d_sf);
    elseif upchirp_peak + downchirp_peak > fft_x_zeropadding*1.5
        cfo_bin = upchirp_peak + downchirp_peak - fft_x_zeropadding*2 - 2;
        cfo = -cfo_bin/2/fft_x_zeropadding * d_bw;
        windows_offset = (downchirp_peak - upchirp_peak) / 2^(11-d_sf);
    else
        cfo_bin = upchirp_peak + downchirp_peak - fft_x_zeropadding - 2;
        cfo = -cfo_bin/2/fft_x_zeropadding * d_bw;
        windows_offset = (fft_x_zeropadding - (upchirp_peak - downchirp_peak)) / 2^(11-d_sf);
    end