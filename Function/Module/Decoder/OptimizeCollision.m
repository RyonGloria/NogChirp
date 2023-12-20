% 测试CHchirpDecoder
clear;
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 读取配置和验证文件
sf = 10;
bw = 125e3;
samplesRate = 2e6;
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
loraSet.channelNum = 4; % 信道数目
loraSet.subchirpNum = 4; % subchirp数目
loraSet.payloadNum = 20; % payload数目
DEBUG = false;
% 从文件中加载bin groundtrurh
true_bin = importdata(strcat('.\Code\Config\bin\SF', string(sf), '.txt'))';
% 读取文件夹下所有采样值文件
fileDir = 'E:\CHchirp\Samples\SF' + string(loraSet.sf) ...
    + '\BW' + string(loraSet.bw/1000) ...
    + '\subchirp' + string(loraSet.subchirpNum) ...
    + '\channel' + string(loraSet.channelNum) ...
    + '\';
fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
% 初始化decoder
CHchirpMultiDecoder = CHchirpMultiDecoder(loraSet, DEBUG);
% 设置循环次数
SignalLength = loraSet.dine*200;  % 整个信号的最大长度
times = 100;
pkgNum = 2;  % 冲突包的数目
% 做times次冲突仿真
for time = 1:times
    % 合成冲突
    fileIndex = randi([1, length(fileIn)], 1, pkgNum);
    winoff = randi([1, 30], 1, pkgNum); % 随机窗口
    off = randi([1, loraSet.dine], 1, pkgNum); % 窗口内随机off
    power = rand(1, pkgNum);
    signalAll = zeros(1, SignalLength);
    % 从文件中随机读取pkgNum个信号,设置随机偏移，合成冲突信号
    for pkgIndex = 1:pkgNum
        signal = readSignalFile(fileDir, fileIn(fileIndex(pkgIndex)));
        offset = loraSet.dine*winoff(pkgIndex) + off(pkgIndex);
        signal = [zeros(1, offset), signal, zeros(1, SignalLength - length(signal) - offset)];
        signalAll = signalAll + signal*power(pkgIndex);
    end
    % 解码
    CHchirpMultiDecoder = CHchirpMultiDecoder.decode(signalAll);
    % 验证结果
    trueAll = 0;
    for binResultIndex = 1:length(CHchirpMultiDecoder.binArray)
        [true_chirp, true_rate] = vertify_bin(CHchirpMultiDecoder.binArray{binResultIndex}, true_bin);
        if true_rate == 1
            trueAll = trueAll + 1/pkgNum;
            if trueAll == 1
                break;
            end
        end
    end
    disp("冲突实验序号:" + time + ", 正确率为：" + trueAll*100 + "%");
end

toc;
fclose all;