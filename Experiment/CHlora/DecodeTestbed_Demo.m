% 画出所选文件的STFT，并给出bin值
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 采样值文件读取路径和保存路径
inDir = 'E:\share\samples\';
% 读取配置和验证文件
[loraSet] = readLoraSet('sf10_BW125.json', 2e6);
[loraSetSF6] = readLoraSet('sf6_BW31_25.json', 2e6);
bw = loraSet.bw;
gbw = bw/2;
preambleChannelChoice = [-(bw+gbw)/2-(bw+gbw), -(bw+gbw)/2, (bw+gbw)/2, (bw+gbw)/2+(bw+gbw)];
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
    % align subchirp
    subchirp = signal(15.25*loraSet.dine+1:15.5*loraSet.dine);
%     stft(subchirp, loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
    [d_downchirp_cfo_SF6, d_upchirp_cfo_SF6] = rebuild_idealchirp_cfo(loraSetSF6, -cfo, 0);
    dechirp = subchirp .* d_upchirp_cfo_SF6;
    dechirp_fft = abs(fft(dechirp, loraSetSF6.dine));
    dechirp_fft = dechirp_fft(1:loraSetSF6.fft_x) + dechirp_fft(loraSetSF6.dine-loraSetSF6.fft_x+1:loraSetSF6.dine);
    [~, downchirp_bin] = max(dechirp_fft);

%     cor = abs(xcorr(subchirp, d_downchirp_cfo_SF6, 1000));
%     plot(cor);

%     subchirp = signal(15.25*loraSet.dine+1:15.5*loraSet.dine);
%     stft(subchirp, loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
%     record = zeros(1, loraSetSF6.dine);
%     record_bin = zeros(1, loraSetSF6.dine);
%     for downchirp_samples = 1:loraSetSF6.factor:loraSetSF6.dine
% %         upchirp_tmp = [d_upchirp_cfo_SF6(1:downchirp_samples), zeros(1, loraSetSF6.dine-downchirp_samples)];
% %         upchirp_tmp = [d_downchirp_cfo_SF6(1:downchirp_samples), zeros(1, loraSetSF6.dine-downchirp_samples)];
%         chirp_tmp = [subchirp(1:downchirp_samples), zeros(1, loraSetSF6.dine-downchirp_samples)];
%         cor = abs(xcorr(chirp_tmp, d_downchirp_cfo_SF6, 0));
%         record(downchirp_samples) = cor;
% %         dechirp = chirp_tmp .* d_upchirp_cfo_SF6;
% %         dechirp_fft = abs(fft(dechirp, loraSetSF6.dine));
% %         dechirp_fft = dechirp_fft(1:loraSetSF6.fft_x) + dechirp_fft(loraSetSF6.dine-loraSetSF6.fft_x+1:loraSetSF6.dine);
% %         [record(downchirp_samples), record_bin(downchirp_samples)] = max(dechirp_fft);
%     end
%     plot(record);
%     xx = 1;

    samples_shift = loraSetSF6.dine - downchirp_bin * loraSetSF6.factor; 
    subchirp = signal(15.25*loraSet.dine+1-samples_shift:15.5*loraSet.dine-samples_shift);
%     stft(subchirp, loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
    dechirp = subchirp .* d_upchirp_cfo_SF6;
    dechirp_fft = abs(fft(dechirp, loraSetSF6.dine));
    dechirp_fft = dechirp_fft(1:loraSetSF6.fft_x) + dechirp_fft(loraSetSF6.dine-loraSetSF6.fft_x+1:loraSetSF6.dine);
    [~, downchirp_bin] = max(dechirp_fft);
    plot(dechirp_fft);

%     subchirp = signal(15.25*loraSet.dine+1:15.5*loraSet.dine);
%     stft(subchirp, loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
    record = zeros(1, loraSetSF6.dine);
    record_bin = zeros(1, loraSetSF6.dine);
    for downchirp_samples = 1:loraSetSF6.dine
%         upchirp_tmp = [d_upchirp_cfo_SF6(1:downchirp_samples), zeros(1, loraSetSF6.dine-downchirp_samples)];
%         upchirp_tmp = [d_downchirp_cfo_SF6(1:downchirp_samples), zeros(1, loraSetSF6.dine-downchirp_samples)];
        chirp_tmp = [subchirp(1:downchirp_samples), zeros(1, loraSetSF6.dine-downchirp_samples)];
%         cor = abs(xcorr(chirp_tmp, d_downchirp_cfo_SF6, 0));
%         record(downchirp_samples) = cor;
        dechirp = chirp_tmp .* d_upchirp_cfo_SF6;
        dechirp_fft = abs(fft(dechirp, loraSetSF6.dine));
        dechirp_fft = dechirp_fft(1:loraSetSF6.fft_x) + dechirp_fft(loraSetSF6.dine-loraSetSF6.fft_x+1:loraSetSF6.dine);
        [record(downchirp_samples), record_bin(downchirp_samples)] = max(dechirp_fft);
    end
    plot(record);
    xx = 1;

% 
%     % search for max peak
%     record_bin = zeros(0);
%     for samShift = -loraSetSF6.factor*2 : loraSetSF6.factor*2
%         SFD_signal = signal(15.25*dine+1 - samples_shift + samShift: 15.5*dine - samples_shift + samShift);
%         dechirp = SFD_signal .* d_upchirp_cfo_SF6;
%         dechirp_fft = abs(fft(dechirp, loraSetSF6.dine));
%         dechirp_fft = dechirp_fft(1:loraSetSF6.fft_x) + dechirp_fft(loraSetSF6.dine-loraSetSF6.fft_x+1:loraSetSF6.dine);
%         [max_peak, max_bin] = max(dechirp_fft);
%         if max_bin == 1
%             record_bin = [record_bin; [samShift, max_peak]];
%         end
%     end
%     [~, index] = max(record_bin(:, 2));
%     samShift = record_bin(index, 1);
%     samples_shift_tmp = samples_shift - samShift;
% 
%     subchirp = signal(15.25*loraSet.dine+1-samples_shift_tmp:15.5*loraSet.dine-samples_shift_tmp);
% %     stft(subchirp, loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
%     dechirp = subchirp .* d_upchirp_cfo_SF6;
%     dechirp_fft = abs(fft(dechirp, loraSetSF6.dine));
%     dechirp_fft = dechirp_fft(1:loraSetSF6.fft_x) + dechirp_fft(loraSetSF6.dine-loraSetSF6.fft_x+1:loraSetSF6.dine);
% %     plot(dechirp_fft);
%     [~, downchirp_bin] = max(dechirp_fft);
% 
%     subchirp = signal(15.25*loraSet.dine+1-samples_shift_tmp:15.5*loraSet.dine-samples_shift_tmp);
%     stft(subchirp, loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
%     record = zeros(1, loraSetSF6.dine);
%     record_bin = zeros(1, loraSetSF6.dine);
%     for downchirp_samples = 1:loraSetSF6.dine
%         upchirp_tmp = [d_upchirp_cfo_SF6(1:downchirp_samples), zeros(1, loraSetSF6.dine-downchirp_samples)];
%         dechirp = subchirp .* upchirp_tmp;
%         dechirp_fft = abs(fft(dechirp, loraSetSF6.dine));
%         dechirp_fft = dechirp_fft(1:loraSetSF6.fft_x) + dechirp_fft(loraSetSF6.dine-loraSetSF6.fft_x+1:loraSetSF6.dine);
%         [record(downchirp_samples), record_bin(downchirp_samples)] = max(dechirp_fft);
%     end
%     plot(record);

    signal_tmp = singalOut(:, 15.25*loraSet.dine+3329+1:end);
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