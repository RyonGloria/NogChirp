function [channelPredictArray] = predictChannelVertex(loraSet, signal, downchirpCfo, chirpOffset, chirpNum, subchirpNum)
    fftX = loraSet.fft_x;
    dine = loraSet.dine;
    channelPredictArray = zeros(chirpNum, subchirpNum);
    for chirpCount = 0:chirpNum - 1
        subchirpChannel = zeros(1,subchirpNum);
        for subchirpCount = 1:subchirpNum
            channelTmp = zeros(1,subchirpNum);
            for channelPredict = 1:4
                chirpIntegrated = signal(channelPredict, (chirpOffset+chirpCount)*dine+(subchirpCount-1)*dine/subchirpNum+1 : (chirpOffset+chirpCount)*dine+subchirpCount*dine/subchirpNum);
                dechirp = downchirpCfo((subchirpCount-1)*dine/subchirpNum+1 : subchirpCount*dine/subchirpNum) .* chirpIntegrated;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
                channelTmp(channelPredict) = max(dechirp_fft);
            end
            [~, subchirpChannel(subchirpCount)] = max(channelTmp);
        end
        channelPredictArray(chirpCount+1, :) = subchirpChannel;
    end