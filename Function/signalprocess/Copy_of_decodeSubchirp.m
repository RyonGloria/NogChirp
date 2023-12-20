function [chirp1Bin] = Copy_of_decodeSubchirp(loraSet, signal, channelChoice, channelChoiceIndexTable, downchirpCfo, chirpOffset, chirpNum, subchirpNum, channelChoiceNum, upchirpCfo)
    fftX = loraSet.fft_x;
    dine = loraSet.dine;
    channelChoice = mod(channelChoice, channelChoiceNum);
    chirp1Bin = zeros(1, chirpNum);
    for chirpCount = 0:chirpNum-1
%         dechirp_fft_merge = zeros(1,fftX);
        dechirp_bin_tmp = zeros(1, subchirpNum);
        for subchirpCount = 0:subchirpNum-1
            subplot(4,4,subchirpCount+1);
            channelChoiceTmp = channelChoice(chirpCount*subchirpNum+1+subchirpCount)+1;
            chirpIntegrated = signal(channelChoiceIndexTable(channelChoiceTmp), (chirpOffset+chirpCount)*dine+subchirpCount*dine/subchirpNum+1 : (chirpOffset+chirpCount)*dine+(subchirpCount+1)*dine/subchirpNum);
%             if subchirpCount <= 7
%                 chirpIntegrated = chirpIntegrated + upchirpCfo(1:dine/subchirpNum);
%             end
            dechirp = downchirpCfo(subchirpCount*dine/subchirpNum+1 : (subchirpCount+1)*dine/subchirpNum) .* chirpIntegrated;
%             dechirp_fft = abs(fft(dechirp, dine));
%             dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(end-fftX+1:end);
%             max_peak = 0;
%             for pahse = 0:2*pi/180*128:2*pi
%                 dechirp_fft_result = dechirp_fft(1:fftX).*exp(1i*pahse) + dechirp_fft(dine-fftX+1:dine);
%                 dechirp_fft_result_abs = abs(dechirp_fft_result);
%                 max_value = max(dechirp_fft_result_abs);
%                 if max_value > max_peak
%                     max_peak = max_value;
%                     dechirp_fft_tmp = dechirp_fft_result_abs;
%                 end
%             end
%             dechirp_fft = dechirp_fft_tmp;
            [dechirp_fft] = phaseCompensate(dechirp, loraSet);
% plot(dechirp_fft);

            [~, dechirp_bin_tmp(subchirpCount+1)] = max(dechirp_fft);
            plot(dechirp_fft);
%             stft(dechirp, loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
            [amp, pos] = max(dechirp_fft);
            fprintf("pos: %f, peak: %f\n", pos, amp);
        end
%         [~, chirp1Bin(chirpCount+1)] = max(dechirp_fft_merge);
        chirp1Bin(chirpCount+1) = mode(dechirp_bin_tmp);
    end


        % for subNum = 1:4
%     subplot(2,2,subNum);
%     dechirp = samples(dine*(subNum-1)/4+1:dine*subNum/4) .* downchirp(dine*(subNum-1)/4+1:dine*subNum/4);
%     dechirp_fft = fft(dechirp, dine);
%     dechirp_fft_tmp = [];
%     max_peak = 0;
%     for pahse = 0:2*pi/180:2*pi
%         dechirp_fft_result = dechirp_fft(1:fftX).*exp(1i*pahse) + dechirp_fft(dine-fftX+1:dine);
%         max_value = max(abs(dechirp_fft_result));
%         if max_value > max_peak
%             max_peak = max_value;
%             dechirp_fft_tmp = dechirp_fft_result;
%         end
%     end
%     plot(abs(dechirp_fft_tmp));
% end