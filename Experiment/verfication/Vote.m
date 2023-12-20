% 测试投票机制
clear;
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 采样值文件读取路径和保存路径
% inDir = 'E:\share\samples\';
% 读取配置和验证文件
sf = 9;
bw = 125e3;


samplesRate = 2e6;
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
loraSet.channelNum = 2; % 信道数目
loraSet.subchirpNum = 2; % subchirp数目
loraSet.payloadNum = 20; % payload数目
% 从文件中加载bin groundtrurh
true_bin = importdata(strcat('.\Code\Config\bin\SF', string(sf), '.txt'))';
% 读取文件夹下所有采样值文件
fileDir = 'E:\CHchirp\Samples\SF' + string(loraSet.sf) ...
    + '\BW' + string(loraSet.bw/1000) ...
    + '\subchirp' + string(loraSet.subchirpNum) ...
    + '\channel' + string(loraSet.channelNum) ...
    + '\';
fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
CHchirpDecoder = CHchirpDecoder(loraSet);
% [signal2] = readSignalFile(inDir, fileIn(2));
% [signal8] = readSignalFile(inDir, fileIn(8));
SignalLength = loraSet.dine*200;
% 添加offset，合成信号
for off = 1:1000
    winoff = randi([1, 10]); % 随机窗口
    bigOff = randi([1, 11]);
    % 从文件中随机读取两个信号
    file1Index = randi([1, length(fileIn)]);
    file2Index = randi([1, length(fileIn)]);
    [signal2] = readSignalFile(fileDir, fileIn(file1Index));
    [signal8] = readSignalFile(fileDir, fileIn(file2Index));
    offset = loraSet.dine*winoff + fix(loraSet.dine * bigOff/11) + off;
    signal2 = [signal2, zeros(1, SignalLength-length(signal2))];
    signal8 = [zeros(1, offset), signal8, zeros(1, SignalLength - length(signal8) - offset)];
    signal_all = signal2 + signal8*0.5;
    % 解调信号
    CHchirpDecoder = CHchirpDecoder.decode(signal_all);
    % 验证解调出来的所有bin是否正确
    [true_chirp, true_rate] = vertify_bin(CHchirpDecoder.payloadBin, true_bin);
    if true_rate ~= 1
        disp("[参数]: " + "SF: " + string(loraSet.sf) ...
            + "channelNum: " + string(loraSet.channelNum) ...
            + "subchirpNum: " + string(loraSet.subchirpNum));
        disp("[decode]: "+" off为：" + off + "正确率为：" + true_rate*100 + "%");
        break;
    end
end

CHchirpDecoder = CHchirpDecoder.decodeVote(signal_all);
% 验证解调出来的所有bin是否正确
[true_chirp, true_rate] = vertify_bin(CHchirpDecoder.payloadBin, true_bin);
disp("[decodeVote]: "+"off为：" + off + "正确率为：" + true_rate*100 + "%");

toc;
fclose all;