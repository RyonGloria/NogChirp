% 测试CHchirpDecoder
clear;
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 读取配置和验证文件
sf = 10;
bw = 250e3;
samplesRate = 2e6;
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
loraSet.channelNum = 6; % 信道数目
loraSet.subchirpNum = 2; % subchirp数目
loraSet.payloadNum = 20; % payload数目
DEBUG = false;
% 从文件中加载bin groundtrurh
true_bin = importdata(strcat('.\Code\Config\bin\SF', string(sf), '.txt'))';
% 读取文件夹下所有采样值文件
fileDir = 'D:\CHchirp23_11\SF' + string(loraSet.sf) ...
        + '_BW' + string(loraSet.bw/1000) ...
        + 'K_sub' + string(loraSet.subchirpNum) ...
        + '_channels' + string(loraSet.channelNum) ...
        + '\';
fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
CHchirpMultiDecoder = CHchirpMultiDecoder(loraSet, DEBUG);
record = cell(3, 1);  % 设置元胞数组记录结果
for fileCount = 1:length(fileIn)
% for fileCount = 59:59
    % 从文件中读取信号流
    [signal] = readSignalFile(fileDir, fileIn(fileCount));

%     figure(1);
%     stft(signal(1:40*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);

    CHchirpMultiDecoder = CHchirpMultiDecoder.decode(signal);
    
    % 验证解调出来的所有bin是否正确
    trueRateMax = 0;
    for binResultIndex = 1:length(CHchirpMultiDecoder.binArray)
        [true_chirp, true_rate] = vertify_bin(CHchirpMultiDecoder.binArray{binResultIndex}, true_bin);
        trueRateMax = max(trueRateMax, true_rate);
    end
    if DEBUG
        disp("decode成功的数目：" + length(CHchirpMultiDecoder.binArray)); 
    end
    
    disp("文件序号:" + fileCount + ", 正确率为：" + trueRateMax*100 + "%");
%     record(fileCount) = true_rate;
end

% binEdges = 0:0.01:1;  % 这里设置了20个柱子，可以根据需要调整
% histogram(record, binEdges, 'Normalization', 'probability');
% title('成功率分布');
% xlabel('成功率');
% ylabel('概率');


toc;
fclose all;