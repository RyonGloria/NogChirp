% 测试CHchirpDecoder
clear;
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 读取配置和验证文件
sf = 10;
bw = 250e3;
samplesRate = 2e6;
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
loraSet.channelNum = 4; % 信道数目
loraSet.subchirpNum = 4; % subchirp数目
loraSet.payloadNum = 20; % payload数目
% 从文件中加载bin groundtrurh
true_bin = importdata(strcat('.\Code\Config\bin\SF', string(sf), '.txt'))';
% 读取文件夹下所有采样值文件
fileDir = 'D:\CHchirp23_11\SF' + string(loraSet.sf) ...
        + '_BW' + string(loraSet.bw/1000) ...
        + 'K_sub' + string(loraSet.subchirpNum) ...
        + '_channels' + string(loraSet.channelNum) ...
        + '\';
fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
CHchirpDecoder = CHchirpDecoder(loraSet);

for fileCount = 1:length(fileIn)
% for fileCount = 24:24
    % 从文件中读取信号流
    [signal] = readSignalFile(fileDir, fileIn(fileCount));
%     figure(1);
%     stft(signal(1:40*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);

    try
        CHchirpDecoder = CHchirpDecoder.decode(signal);
    catch
        CHchirpDecoder.payloadBin = zeros(1, loraSet.payloadNum);
    end
    
    % 验证解调出来的所有bin是否正确
    [true_chirp, true_rate] = vertify_bin(CHchirpDecoder.payloadBin, true_bin);
    disp("文件序号:" + fileCount + ", 正确率为：" + true_rate*100 + "%");
    % 如果正确率不为100%，则删除信号
    if true_rate ~= 1
        filePath = strcat(fileDir, fileIn(fileCount).name);
        delete(filePath);
    end
end

% binEdges = 0:0.01:1;  % 这里设置了20个柱子，可以根据需要调整
% histogram(record, binEdges, 'Normalization', 'probability');
% title('成功率分布');
% xlabel('成功率');
% ylabel('概率');


toc;
fclose all;