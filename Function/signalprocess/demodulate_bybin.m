function [subchirp_bin] = demodulate_bybin(loraSet, signal, channel_jump_array, subchirpNum, d_downchirp_cfo)
    subchirp_bin = zeros(1, subchirpNum);
    fftX = loraSet.fft_x;
    dine = loraSet.dine;
    for subchirpCount = 0:subchirpNum-1
        chirpIntegrated = signal(channel_jump_array(subchirpCount+1), subchirpCount*dine/subchirpNum+1 : (subchirpCount+1)*dine/subchirpNum);
        dechirp = d_downchirp_cfo(subchirpCount*dine/subchirpNum+1 : (subchirpCount+1)*dine/subchirpNum) .* chirpIntegrated;
        dechirp_fft = abs(fft(dechirp, dine));
        dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
        [~, subchirp_bin(subchirpCount+1)] = max(dechirp_fft);
    end