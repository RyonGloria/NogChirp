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
% channel_jump_array = zeros(16,16*4);
% channel1 = 1; channel2 = 1; channel3 = 1; channel4 = 0;
% for i = 1:16
%     for k = 0:15
%         channel4 = channel4 + 1;
%         if channel4 == 5
%             channel3 = channel3 + 1;
%             channel4 = 1;
%             if channel3 == 5
%                 channel2 = channel2 + 1;
%                 channel3 = 1;
%                 if channel2 == 5
%                     channel1 = channel1 + 1;
%                     channel2 = 1;
%                 end
%             end
%         end
%         channel_jump_array(i, k*4+1) = channel1;
%         channel_jump_array(i, k*4+2) = channel2;
%         channel_jump_array(i, k*4+3) = channel3;
%         channel_jump_array(i, k*4+4) = channel4;
%     end
% end
% channel_jump_array = [2,1,3,3];
% channel_jump_array = [1,4,3,2,2,1,3,3,3,1,2,2,3,2,3,2,2,1,2,4,4,3,1,1,4,2,2,1,4,4,3,3,3,3,2,3,4,3,4,1,2,4,1,3,1,1,4,2,4,1,3,3,2,4,1,2,4,3,3,4,3,1,1,2,4,4,3,2,2,1,1,2];
% 生成preamble对应信道的idealchirp
[downchirp, upchirp] = buildIdealchirp(loraSet, 0); % build idealchirp
% 读取文件夹下所有采样值文件
fileIn = dir(fullfile(inDir, '*.sigmf-data'));          
for fileCount = 1:length(fileIn)
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

    
    signal_tmp = singalOut(:, 15.25*loraSet.dine-620+1:end);
    subchirpNum = 4;
    for chirp_count = 0:15
        dechirp_bin_tmp = zeros(1, subchirpNum);
        for subchirpCount = 0:subchirpNum-1
%             chirpIntegrated = signal_tmp(channel_jump_array(chirp_count*4+subchirpCount+1), chirp_count*dine + subchirpCount*dine/subchirpNum+1 : chirp_count*dine + (subchirpCount+1)*dine/subchirpNum);
%             chirpIntegrated = signal_tmp(channel_jump_array(fileCount, chirp_count*4+subchirpCount+1), chirp_count*dine + subchirpCount*dine/subchirpNum+1 : chirp_count*dine + (subchirpCount+1)*dine/subchirpNum);
            chirpIntegrated = signal_tmp(channel_jump_array(chirp_count*4+subchirpCount+1), chirp_count*dine + subchirpCount*dine/subchirpNum+1 : chirp_count*dine + (subchirpCount+1)*dine/subchirpNum);
            dechirp = d_downchirp_cfo(subchirpCount*dine/subchirpNum+1 : (subchirpCount+1)*dine/subchirpNum) .* chirpIntegrated;
            dechirp_fft = abs(fft(dechirp, dine));
            dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
%             subplot(2,2,subchirpCount+1);
%             plot(dechirp_fft);
            [~, dechirp_bin_tmp(subchirpCount+1)] = max(dechirp_fft);
        end
%         mode_tmp = mode(dechirp_bin_tmp);
%         sum_tmp = sum(dechirp_bin_tmp == mode_tmp);
%         if sum_tmp ~= 4
%             fprintf("%d ", channel_jump_array(fileCount, chirp_count*4+1 : chirp_count*4+4));
%             fprintf("\n");
%         end
        fprintf("%d ", dechirp_bin_tmp);
        fprintf("\n");
    end

end

toc;
fclose all;