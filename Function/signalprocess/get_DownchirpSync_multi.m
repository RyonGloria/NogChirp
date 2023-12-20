function [channel_jump_array] = get_DownchirpSync_multi(loraSet, signal, d_upchirp_cfo, amp_ref)

    fftX = loraSet.fft_x;
    dine = loraSet.dine;
    dechirp = signal .* d_upchirp_cfo;
    dechirp_fft = abs(fft(dechirp, dine));
    dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
    [pos_array, peak_array] = fft_filterN_peak(loraSet, dechirp_fft, 40);
    dif_array = abs(peak_array - amp_ref);
    [~, min_pos] = min(dif_array);
    downchirp_bin = pos_array(min_pos);
%     [tmp, downchirp_bin] = max(dechirp_fft);
    load("C:\Users\ZKevin\Desktop\CHchirp\random_record.mat");
    channel_jump_array = random_record(downchirp_bin, :);