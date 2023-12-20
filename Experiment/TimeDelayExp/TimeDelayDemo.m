% 画出所选文件的STFT，并给出bin值
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 采样值文件读取路径和保存路径
inDir = 'E:\CHlora_samples\TimeDelay\';
% 读取配置和验证文件
[loraSet] = readLoraSet('sf10_BW125.json', 2e6);
verBin = load("config\TimeDelay.txt");
subchirpNum = 4;
chirp_num = 16;
bw = loraSet.bw;
gbw = bw/2;
preambleChannelChoice = [-(bw+gbw)/2-(bw+gbw), -(bw+gbw)/2, (bw+gbw)/2, (bw+gbw)/2+(bw+gbw)];
% 生成preamble对应信道的idealchirp
[downchirp, upchirp] = buildIdealchirp(loraSet, 0); % build idealchirp
% 读取文件夹下所有采样值文件
fileIn = dir(fullfile(inDir, '*.sigmf-data'));    
result_array = zeros(1, 401);
for fileCount = 1:length(fileIn)
    if mod(fileCount, 10) == 0           % 每循环100次，输出当前进度
        fprintf("The time is %d\n", fileCount);
        toc;
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
    record_shift = zeros(0);
    for samples_shift = -200:200
        signal_tmp = singalOut(:, 15.25*loraSet.dine-620+1 + samples_shift : end);
        true_flag = 1;
        for chirp_count = 0:15
            if true_flag == 0   % when one chirp bin is error, break
                break;
            end
            dechirp_bin_tmp = zeros(1, subchirpNum);
            dechirp_value_tmp = zeros(1, subchirpNum);
            for subchirpCount = 0:subchirpNum-1
                chirpIntegrated = signal_tmp(channel_jump_array(chirp_count*4+subchirpCount+1), chirp_count*dine + subchirpCount*dine/subchirpNum+1 : chirp_count*dine + (subchirpCount+1)*dine/subchirpNum);
                dechirp = d_downchirp_cfo(subchirpCount*dine/subchirpNum+1 : (subchirpCount+1)*dine/subchirpNum) .* chirpIntegrated;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
                [dechirp_value_tmp(subchirpCount+1), dechirp_bin_tmp(subchirpCount+1)] = max(dechirp_fft);
            end
            if sum(dechirp_bin_tmp == verBin(chirp_count+1)) == subchirpNum
                true_flag = true_flag + 1;
            else
                true_flag = 0;
            end
%             fprintf("%d ", dechirp_bin_tmp);
%             fprintf("\n");
        end
        if true_flag == chirp_num + 1
            record_shift = [record_shift; [samples_shift, dechirp_value_tmp(end)]];
        end
    end
    if ~isempty(record_shift)
        [~, time_shift_record_index] = max(record_shift(:, 2));
        time_shift_record = record_shift(time_shift_record_index, 1);
        result_array(time_shift_record+201) = result_array(time_shift_record+201) + 1;
    end
end
plot(-200:200,result_array);
save('TimeDelay_t2000.mat','result_array');

toc;
fclose all;