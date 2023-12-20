% 画出所选文件的STFT，并给出bin值
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 采样值文件读取路径和保存路径
inDir = 'E:\share\samples\';
% 读取配置和验证文件
[loraSet] = readLoraSet('sf10_BW125.json', 2e6);

% 生成preamble对应信道的idealchirp
% [downchirp, upchirp] = buildIdealchirp(loraSet, 375e3); % build idealchirp
[downchirp, upchirp] = buildIdealchirp(loraSet, 0); % build idealchirp
% 读取文件夹下所有采样值文件
fileIn = dir(fullfile(inDir, '*.sigmf-data'));          
for fileCount = 1:1
    if mod(fileCount, 100) == 0           % 每循环100次，输出当前进度
        fprintf("The time is %d\n", fileCount);
    end
    [signal] = readSignalFile(inDir, fileIn(fileCount));
    % find Preamble
%     [PreambleStartPos] = detect_preamble(signal, loraSet, downchirp); 
%     if PreambleStartPos == 999
%         continue;
%     elseif PreambleStartPos ~= 1
%         signal = circshift(signal, -(PreambleStartPos-1) * loraSet.dine);
%     end
%     % 调整cfo和timeoffset
%     [cfo, windowsOffset] = get_cfo_winoff(signal, loraSet, downchirp, upchirp, loraSet.factor, true);
%     signal = circshift(signal, -round(windowsOffset));
    stft(signal(1:60*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
%     stft(signal(15.25*loraSet.dine+1:16.25*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);

%     [downchirp_3175_cfo, upchirp_3175_cfo] = rebuild_idealchirp_cfo(loraSet, cfo, 375e3);
%     tmp = signal(1:loraSet.dine) .* downchirp_3175_cfo;
%     fft_x = loraSet.fft_x;
%     dine = loraSet.dine;
%     dechirp_1_fft = abs(fft(tmp));
%     dechirp_1_fft = dechirp_1_fft(1:fft_x) + dechirp_1_fft(dine-fft_x+1:dine);
%     [~, chirp1_bin] = max(dechirp_1_fft);
%     disp(chirp1_bin);
end

toc;
fclose all;