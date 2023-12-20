function [d_downchirp, d_upchirp] = buildIdealchirp(lora_set, f_temp)
    cmx = 1+1*1i;
    pre_dir = 2*pi;
    f0 = f_temp+lora_set.bw/2;                           % 设置理想upchirp和downchirp的初始频率
    f1 = -f_temp+lora_set.bw/2; 
    d_symbols_per_second = lora_set.bw / lora_set.fft_x; 
    T = -0.5 * lora_set.bw * d_symbols_per_second;  
    d_samples_per_second = lora_set.sample_rate;        % sdr-rtl的采样率
    d_dt = 1/d_samples_per_second;         % 采样点间间隔的时间
    t = d_dt*(0:1:lora_set.dine-1);
    % 计算理想downchirp和upchirp存入d_downchirp和d_upchirp数组中（复数形式）
    d_downchirp = cmx * (cos(pre_dir .* t .* (f0 + T * t)) + sin(pre_dir .* t .* (f0 + T * t))*1i);
    d_upchirp = cmx * (cos(pre_dir .* t .* (f1 + T * t) * -1) + sin(pre_dir .* t .* (f1 + T * t) * -1)*1i);
%     d_upchirp = cmx * (cos(pre_dir .* t .* (-f1 + T * t) * -1) + sin(pre_dir .* t .* (-f1 + T * t) * -1)*1i);