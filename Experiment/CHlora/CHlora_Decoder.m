% 画出所选文件的STFT，并给出bin值
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 采样值文件读取路径和保存路径
inDir = 'E:\share\samples\';
% 读取配置和验证文件
[loraSet] = readLoraSet('sf10_BW125.json', 2e6);
bw = loraSet.bw;  % preamble占用的带宽
gbw = bw/2;     % 过渡带带宽
subchirpNum = 4;    % subchirp的数目
% 每个信道对应的中心频率（目前的信号已经进过频移至0频，其中心频率为采样时设置的中心频率）
preambleChannelChoice = [-(bw+gbw)/2-(bw+gbw), -(bw+gbw)/2, (bw+gbw)/2, (bw+gbw)/2+(bw+gbw)];
% 生成中心频率为0的idealchirp
[downchirp, upchirp] = buildIdealchirp(loraSet, 0); % build idealchirp
% 信号bin的groundtruth
true_bin = [0, 56, 112, 168, 224, 280, 336, 392, 448, 504, 560, 616, 672, 728, 784, 840] + 1;
% 读取文件夹下所有采样值文件
fileIn = dir(fullfile(inDir, '*.sigmf-data'));          
% for fileCount = 1:length(fileIn)
for fileCount = 1:1
    % 从文件中读取信号流
    [signal] = readSignalFile(inDir, fileIn(fileCount));
    % find active channel
    % 将信号拆分成4个信道信号
    [signalOut] = divideChannel(loraSet, signal, preambleChannelChoice, false);
    % 检测4个信道中最有可能存在preamble的信道（FFT幅值最大的）
    [max_amp, max_bin, active_channel] = detect_active_channel(signalOut, loraSet, downchirp);
    % 获得确定信道的信号
    signal_detect = signalOut(active_channel, :);   
    % 检测preamble，确定存在preamble并且获得第一个preamble出现的窗口
    [PreambleStartPos] = detect_preamble(signal_detect, loraSet, downchirp); 
    if PreambleStartPos == 999
        quit;  % 找不到preamble，直接退出程序
    elseif PreambleStartPos ~= 1  % 调整到preamble的第一个窗口
        signal_detect = circshift(signal_detect, -(PreambleStartPos-1) * loraSet.dine);
        signal = circshift(signal, -(PreambleStartPos-1) * loraSet.dine);
    end
    % 通过preamble和SFD的bin来计算CFO和winoffset
    [cfo, windowsOffset] = get_cfo_winoff(signal_detect, loraSet, downchirp, upchirp, loraSet.factor, true);
    % 调整信号的winoffset
    signal = circshift(signal, -round(windowsOffset));
    signal_detect = circshift(signal_detect, -round(windowsOffset));
    % 根据cfo重新生成带有decfo的idealchirp，用于解调
    [d_downchirp_cfo, d_upchirp_cfo] = rebuild_idealchirp_cfo(loraSet, cfo, 0);
    
    % 获得payload阶段的跳信道矩阵
    [channel_jump_array] = getDownchirpSync(loraSet, signal_detect(14.25*loraSet.dine+1:15.25*loraSet.dine), d_upchirp_cfo);
    
    % 将信号四个信道滤波后保存在四个流里面
    [singalOut] = divideChannel(loraSet, signal, preambleChannelChoice, false);

    % 根据第一个跳信道subchirp（bin为0）来重新对齐
    [time_off] = align_windows_bysubchirp(loraSet, singalOut(:, 15.25*loraSet.dine+1:16.25*loraSet.dine), channel_jump_array, subchirpNum, d_downchirp_cfo);

    % demodulate 解调信号
    [subchirp_bin] = demodulate_subchirp(loraSet, singalOut(:, 15.25*loraSet.dine-time_off+1:end), channel_jump_array, subchirpNum, d_downchirp_cfo, 16);
    
    % 验证解调出来的所有bin是否正确
    [true_chirp, true_rate] = vertify_bin(subchirp_bin, true_bin);
    disp("正确率为：" + true_rate*100 + "%");
end

toc;
fclose all;