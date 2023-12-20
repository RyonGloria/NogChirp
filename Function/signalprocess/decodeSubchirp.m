function [chirp1Bin] = decodeSubchirp(loraSet, signal, channelChoice, channelChoiceIndexTable, downchirpCfo, chirpOffset, chirpNum, subchirpNum, channelChoiceNum)
    fftX = loraSet.fft_x;
    dine = loraSet.dine;
    channelChoice = mod(channelChoice, channelChoiceNum);
    chirp1Bin = zeros(1, chirpNum);
    for chirpCount = 0:chirpNum-1
%         dechirp_fft_merge = zeros(1,fftX);
        dechirp_bin_tmp = zeros(1, subchirpNum);
        for subchirpCount = 0:subchirpNum-1
            channelChoiceTmp = channelChoice(chirpCount*subchirpNum+1+subchirpCount)+1;
            chirpIntegrated = signal(channelChoiceIndexTable(channelChoiceTmp), (chirpOffset+chirpCount)*dine+subchirpCount*dine/subchirpNum+1 : (chirpOffset+chirpCount)*dine+(subchirpCount+1)*dine/subchirpNum);
            dechirp = downchirpCfo(subchirpCount*dine/subchirpNum+1 : (subchirpCount+1)*dine/subchirpNum) .* chirpIntegrated;
            dechirp_fft = abs(fft(dechirp, dine));
            dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
            [~, dechirp_bin_tmp(subchirpCount+1)] = max(dechirp_fft);
        end
%         [~, chirp1Bin(chirpCount+1)] = max(dechirp_fft_merge);
        chirp1Bin(chirpCount+1) = mode(dechirp_bin_tmp);
    end