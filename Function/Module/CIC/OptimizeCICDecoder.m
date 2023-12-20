% 测试CHchirpDecoder
clear;
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 读取配置和验证文件
sf = 9;
bw = 250e3;
samplesRate = 2e6;
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
loraSet.channelNum = 4; % 信道数目
loraSet.subchirpNum = 4; % subchirp数目
loraSet.payloadNum = 20; % payload数目
DEBUG = false;
% 从文件中加载bin groundtrurh
true_bin = importdata(strcat('.\Code\Config\bin\SF', string(sf), '.txt'))';
% 读取文件夹下所有采样值文件
fileDir = 'D:\CHchirp_20Nodes\Pyramid_SF9_BW250K\Node1\';
fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
% 初始化decoder
CICDecoder = CICDecoder(loraSet);
for fileCount = 1:length(fileIn)
% for fileCount = 59:59
    % 从文件中读取信号流
    [signal] = readSignalFile(fileDir, fileIn(fileCount));
    signal = signal(1, loraSet.dine*40);
    winoff = randi([2, 30], 1, 1); % 随机窗口
    off = randi([1, loraSet.dine], 1, 1); % 窗口内随机off
    SIR = (rand(1, 1) * 10) - 5;
%     signal = signal + 0.5*circshift(signal, loraSet.dine*20.75);
    CICDecoder = CICDecoder.decode(signal);
    trueRate = calculateAccuracy(true_bin, CICDecoder.binRecord{1});
end