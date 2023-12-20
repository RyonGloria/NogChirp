function [Preamble_start_pos, Preamble_num, Preamble_bin] = detect_preamble_bin(lora_set, signal, downchirp)
    dine = lora_set.dine;
    fft_x = lora_set.fft_x;
    preamble_len = lora_set.Preamble_length;
    candidate = zeros(1, preamble_len + 2);
%     detect_pos = loraSet.fft_x + 2;
    for t = 1:preamble_len + 2
        signal_tmp = signal((t-1)*dine+1 : t*dine);
        dechirp = signal_tmp .* downchirp;
        dechirp_fft = abs(fft(dechirp, dine));
        dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine);
        [~, candidate(t)] = max(dechirp_fft);
    end
    Preamble_bin = mode(candidate);
    Preamble_start_pos = find(candidate == Preamble_bin);
    Preamble_start_pos = Preamble_start_pos(1);
    Preamble_num = 0;
    for t = Preamble_start_pos:preamble_len + 2
        if candidate(t) == Preamble_bin
            Preamble_num = Preamble_num + 1;
        else
            break;
        end
    end