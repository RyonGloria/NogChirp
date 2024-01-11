% 测试CHchirpDecoder
clear;
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 读取配置和验证文件
sf = 10;
bw = 125e3;
samplesRate = 2e6;
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
loraSet.payloadNum = 23; % payload数目
NogChirpDecoder = NogChirpDecoder(loraSet);
% 读取文件夹下所有采样值文件
fileDir = '\\192.168.3.102\e\data\delay_231219\';
fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
% 从文件中读取信号流
[signal] = readSignalFile(fileDir, fileIn(2));

NogChirpDecoder = NogChirpDecoder.decode(signal);

% disp(NogChirpDecoder.payloadBin);

toc;
fclose all;