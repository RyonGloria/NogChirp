function [samples_fft_merge] = get_zeropadding_fft(signal, lora_set, downchirp)
    zeropadding_size = lora_set.factor;                   % 设置补零的数量，这里的decim表示，补上decim-1倍窗口的零，计算FFT时一共是decim倍的窗口（decim+1）
    d_sf = lora_set.sf;
    dine = lora_set.dine;
    fft_x = lora_set.fft_x;

    dine_zeropadding = dine*zeropadding_size*2^(10-d_sf);
    fft_x_zeropadding = fft_x*zeropadding_size*2^(10-d_sf);

    samples_fft = abs(fft(signal .* downchirp, dine_zeropadding, 2));
    samples_fft_merge = samples_fft(:, 1:fft_x_zeropadding) + samples_fft(:,dine_zeropadding-fft_x_zeropadding+1:dine_zeropadding);
