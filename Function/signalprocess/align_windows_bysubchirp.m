function [time_off] = align_windows_bysubchirp(loraSet, signal, channel_jump_array, subchirpNum, d_downchirp_cfo)
    fftX = loraSet.fft_x;
    dine = loraSet.dine;
    dechirp_bin_tmp = zeros(1, subchirpNum);
    for subchirpCount = 0:subchirpNum-1
        chirpIntegrated = signal(channel_jump_array(subchirpCount+1), subchirpCount*dine/subchirpNum+1 : (subchirpCount+1)*dine/subchirpNum);
        dechirp = d_downchirp_cfo(subchirpCount*dine/subchirpNum+1 : (subchirpCount+1)*dine/subchirpNum) .* chirpIntegrated;
        dechirp_fft = abs(fft(dechirp, dine*loraSet.factor));
        dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
        [~, dechirp_bin_tmp(subchirpCount+1)] = max(dechirp_fft);
    end
    time_off = round(mean(dechirp_bin_tmp));