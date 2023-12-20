function [d_downchirp_cfo, d_upchirp_cfo] = rebuild_idealchirp_cfo(lora_set, cfo, f_temp)
    cmx = 1+1*1i;
    pre_dir = 2*pi;
%     f0 = lora_set.bw/2;                           % 设置理想upchirp和downchirp的初始频率
    d_symbols_per_second = lora_set.bw / lora_set.fft_x; 
    T = -0.5 * lora_set.bw * d_symbols_per_second;  
    d_samples_per_second = lora_set.sample_rate;       % sdr-rtl的采样率
    d_dt = 1/d_samples_per_second;         % 采样点间间隔的时间
    t = d_dt*(0:1:lora_set.dine-1);
%     f0 = lora_set.bw/2+cfo;                      % downchirp和upchirp收到CFO的影响是相反的，所以要重新设置两个的初始频率
%     f1 = lora_set.bw/2-cfo;
    f0 = f_temp+lora_set.bw/2+cfo;                           % 设置理想upchirp和downchirp的初始频率
    f1 = -f_temp+lora_set.bw/2-cfo; 

    % 计算理想downchirp和upchirp存入d_downchirp和d_upchirp数组中（复数形式）
    d_downchirp_cfo = cmx * (cos(pre_dir .* t .* (f0 + T * t)) + sin(pre_dir .* t .* (f0 + T * t))*1i);
    d_upchirp_cfo = cmx * (cos(pre_dir .* t .* (f1 + T * t) * -1) + sin(pre_dir .* t .* (f1 + T * t) * -1)*1i);