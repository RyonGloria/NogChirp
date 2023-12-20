function [max_amp, max_bin, active_channel] = detect_active_channel(signal, lora_set, downchirp)
    channel_num = size(signal, 1);
    dine = lora_set.dine;
    fft_x = lora_set.fft_x;
    max_amp = 0;
    max_bin = 0;
    active_channel = 0;
    for channel = 1:channel_num
        signal_tmp = signal(channel, 1*dine+1:2*dine);
        dechirp = signal_tmp .* downchirp;
        dechirp_fft = abs(fft(dechirp, dine));
        dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine);
        [amp, bin] = max(dechirp_fft);
        if amp > max_amp
            max_amp = amp;
            max_bin = bin;
            active_channel = channel;
        end
%         disp([amp, bin]);
    end