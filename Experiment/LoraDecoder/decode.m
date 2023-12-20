% 测试CHchirpDecoder
clear;
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 读取配置和验证文件
sf = 10;
bw = 125e3;
samplesRate = 2e6;
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
loraSet.payloadNum = 25; % payload数目
NogChirpDecoder = NogChirpDecoder(loraSet);
% 读取文件夹下所有采样值文件
fileDir = '\\192.168.3.102\e\share\samples\';

fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
for fileCount = 1:length(fileIn)
% for fileCount = 59:59
    % 从文件中读取信号流
    [signal] = readSignalFile(fileDir, fileIn(fileCount));
%     signal = [zeros(1, 10*loraSet.dine), signal, zeros(1, 10*loraSet.dine)];
    % signal = downsample(signal, 2);
    % 写入读取区
    NogChirpDecoder = NogChirpDecoder.decode(signal);

    disp(NogChirpDecoder.binRecord);
end

toc;
fclose all;