function [active_channel] = get_active_channel_multi(signal, lora_set, downchirp)
    channel_num = size(signal, 1);
    dine = lora_set.dine;
    fft_x = lora_set.fft_x;
    active_channel = [];
    for channel = 1:channel_num
        signal_tmp = signal(channel, 1:dine);
        dechirp = signal_tmp .* downchirp;
        dechirp_fft = abs(fft(dechirp, dine));
        dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine);
        [amp, ~] = max(dechirp_fft);
        if amp > mean(dechirp_fft) * 10
            active_channel = [active_channel, channel];
        end
    end