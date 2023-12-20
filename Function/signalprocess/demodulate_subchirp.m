function [subchirp_bin] = demodulate_subchirp(loraSet, signal, channel_jump_array, subchirpNum, d_downchirp_cfo, chirpNum)
    subchirp_bin = zeros(chirpNum, subchirpNum);
    fftX = loraSet.fft_x;
    dine = loraSet.dine;
    for chirp_count = 0:chirpNum-1
        dechirp_bin_tmp = zeros(1, subchirpNum);
        for subchirpCount = 0:subchirpNum-1
            chirpIntegrated = signal(channel_jump_array(chirp_count*4+subchirpCount+1), chirp_count*dine + subchirpCount*dine/subchirpNum+1 : chirp_count*dine + (subchirpCount+1)*dine/subchirpNum);
%             chirpIntegrated = signal_tmp(array_tmp(subchirpCount+1), chirp_count*dine + subchirpCount*dine/subchirpNum+1 : chirp_count*dine + (subchirpCount+1)*dine/subchirpNum);
            dechirp = d_downchirp_cfo(subchirpCount*dine/subchirpNum+1 : (subchirpCount+1)*dine/subchirpNum) .* chirpIntegrated;
            dechirp_fft = abs(fft(dechirp, dine));
            dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
            [~, dechirp_bin_tmp(subchirpCount+1)] = max(dechirp_fft);
        end
        subchirp_bin(chirp_count+1, :) = dechirp_bin_tmp;
    end