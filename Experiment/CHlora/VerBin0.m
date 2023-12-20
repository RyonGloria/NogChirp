% 画出所选文件的STFT，并给出bin值
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 采样值文件读取路径和保存路径
inDir = 'E:\share\samples\';
% 读取配置和验证文件
[loraSet] = readLoraSet('sf10_BW250.json', 2e6);

% 生成preamble对应信道的idealchirp
% [downchirp, upchirp] = buildIdealchirp(loraSet, 375e3); % build idealchirp
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
    
    dine = loraSet.dine;
    fft_x = loraSet.fft_x;
    fft_windows = 18;
    samples = signal(12.25*loraSet.dine:end);
    samples = reshape(samples(1:fft_windows*dine),[dine, fft_windows]).';
    samples_dechirp = samples .* d_downchirp_cfo;
    samples_fft = abs(fft(samples_dechirp,dine,2));
    samples_fft_merge = [samples_fft(:,1:fft_x/2) + samples_fft(:,dine-fft_x+1:dine-fft_x/2), samples_fft(:,dine-fft_x/2+1:dine) + samples_fft(:,fft_x/2+1:fft_x)];
    [~, max_pos] = max(samples_fft_merge, [], 2);
    [~,F] = mode(max_pos);
    if F == fft_windows
        fprintf("Find true samples: %d", fileIn(fileCount).name);
    end

end

toc;
fclose all;