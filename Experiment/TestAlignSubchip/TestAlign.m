% 画出所选文件的STFT，并给出bin值
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 采样值文件读取路径和保存路径
inDir = 'E:\share\samples\';
% 读取配置和验证文件
[loraSet] = readLoraSet('sf10_BW125.json', 2e6);
bw = loraSet.bw;
gbw = bw/2;
preambleChannelChoice = [-(bw+gbw)/2-(bw+gbw), -(bw+gbw)/2, (bw+gbw)/2, (bw+gbw)/2+(bw+gbw)];
% 生成preamble对应信道的idealchirp
[downchirp, upchirp] = buildIdealchirp(loraSet, 0); % build idealchirp
% 读取文件夹下所有采样值文件
fileIn = dir(fullfile(inDir, '*.sigmf-data'));        
true_bin = [0, 56, 112, 168, 224, 280, 336, 392, 448, 504, 560, 616, 672, 728, 784, 840] + 1;
for fileCount = 1:1
    if mod(fileCount, 100) == 0           % 每循环100次，输出当前进度
        fprintf("The time is %d\n", fileCount);
    end
    [signal] = readSignalFile(inDir, fileIn(fileCount));
    % find Preamble
    [PreambleStartPos] = detect_preamble(signal, loraSet, downchirp); 
    if PreambleStartPos == 999
        continue;
    elseif PreambleStartPos ~= 1
        signal = circshift(signal, -(PreambleStartPos-1) * loraSet.dine);
    end
    % 调整cfo和timeoffset
    [cfo, windowsOffset] = get_cfo_winoff(signal, loraSet, downchirp, upchirp, loraSet.factor, true);
    signal = circshift(signal, -round(windowsOffset));
    [d_downchirp_cfo, d_upchirp_cfo] = rebuild_idealchirp_cfo(loraSet, cfo, 0);
    % get downchirp_sync bin
    fftX = loraSet.fft_x;
    dine = loraSet.dine;
    downchirp_sync_signal = signal(14.25*loraSet.dine+1:15.25*loraSet.dine);
    dechirp = downchirp_sync_signal .* d_upchirp_cfo;
    dechirp_fft = abs(fft(dechirp, dine));
    dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
    [~, downchirp_bin] = max(dechirp_fft);
    load("random_record.mat");
    channel_jump_array = random_record(downchirp_bin, :);
    [singalOut] = divideChannel(loraSet, signal, preambleChannelChoice, false);
    % sf10 - 968   sf11-(-1080)  sf12-(-1520)
%     signal_tmp = singalOut(:, 15.25*loraSet.dine+1:16.25*loraSet.dine);
%     signal_tmp = [signal_tmp(1,1:dine/4), signal_tmp(4,dine/4+1:dine/2), signal_tmp(3,dine/2+1:dine*3/4), signal_tmp(2,dine*3/4+1:dine)];
%     stft(signal_tmp, loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
%     singal_tmp = signal(1:dine);
%     figure(1);
%     stft(signal(15.25*loraSet.dine-620+1:16.25*loraSet.dine-620), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
%     figure(2);
%     stft(signal(16.25*loraSet.dine-620+1:17.25*loraSet.dine-620), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);

    signal_tmp = singalOut(:, 15.25*loraSet.dine-620+1:end);
    subchirpNum = 4;
    samShiftArray = zeros(1,16*4);
    for chirp_count = 0:15
        dechirp_bin_tmp = zeros(1, subchirpNum);
        for subchirpCount = 0:subchirpNum-1
            max_peak = 0;
            for pahse = 0:pi/3600:2*pi
                chirpIntegrated = signal_tmp(channel_jump_array(chirp_count*4+subchirpCount+1), chirp_count*dine + subchirpCount*dine/subchirpNum+1: chirp_count*dine + (subchirpCount+1)*dine/subchirpNum);
                dechirp = d_downchirp_cfo(subchirpCount*dine/subchirpNum+1 : (subchirpCount+1)*dine/subchirpNum) .* chirpIntegrated.*exp(1i*pahse);
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
                [max_value, dechirp_bin_tmp(subchirpCount+1)] = max(dechirp_fft);
%                 if dechirp_bin_tmp(subchirpCount+1) == true_bin(chirp_count+1)
                    if max_value > max_peak
                        max_peak = max_value;
                        record_phase = pahse;
                    end
%                 end
            end
            samShiftArray(chirp_count*4 + subchirpCount + 1) = record_phase;
        end
%         fprintf("%d ", dechirp_bin_tmp);
%         fprintf("\n");
    end
    for chirp_count = 0:15
%     for chirp_count = 2:2
        chirpIntegrated = zeros(0);
        dechirp_bin_tmp = zeros(1, subchirpNum);
        for subchirpCount = 0:3
            chirpIntegrated = [chirpIntegrated, signal_tmp(channel_jump_array(chirp_count*4+subchirpCount+1), chirp_count*dine + subchirpCount*dine/subchirpNum+1 : chirp_count*dine + (subchirpCount+1)*dine/subchirpNum) .* samShiftArray(chirp_count*4 + subchirpCount + 1)];

%             dechirp = d_downchirp_cfo(subchirpCount*dine/subchirpNum+1 : (subchirpCount+1)*dine/subchirpNum) .* chirpIntegrated;
%             dechirp_fft = abs(fft(dechirp, dine));
%             dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
%             [~, dechirp_bin_tmp(subchirpCount+1)] = max(dechirp_fft);
        end
        dechirp = chirpIntegrated .* d_downchirp_cfo(1:dine);
        dechirp_fft = abs(fft(dechirp, dine));
        dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
%         figure(2);
%         plot(dechirp_fft);
        [~, max_bin] = max(dechirp_fft);
        fprintf("%d ", max_bin);
        fprintf("\n");
%         stft(chirpIntegrated, loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
%         fprintf("%d ", dechirp_bin_tmp);
%         fprintf("\n");
    end
%     figure(1);
%     stft(singalOut(1, 1:40*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
%     figure(2);
%     stft(singalOut(2, 1:40*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
%     figure(3);
%     stft(singalOut(3, 1:40*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
%     figure(4);
%     stft(singalOut(4, 1:40*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
%     signal_tmp(1, 1:dine/4)

end

toc;
fclose all;