% 测试CHchirpDecoder
clear;
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 采样值文件读取路径和保存路径
% inDir = 'E:\share\samples\';
% 读取配置和验证文件
sf = 10;
bw = 125e3;
samplesRate = 2e6;
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
loraSet.channelNum = 4; % 信道数目
loraSet.subchirpNum = 4; % subchirp数目
loraSet.payloadNum = 20; % payload数目