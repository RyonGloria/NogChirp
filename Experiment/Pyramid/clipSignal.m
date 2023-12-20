% 测试CHchirpDecoder
clear;
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 读取配置和验证文件
sf = 10;
bw = 125e3;
samplesRate = 2e6;
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
loraSet.channelNum = 2; % 信道数目
loraSet.subchirpNum = 2; % subchirp数目
loraSet.payloadNum = 20; % payload数目
CHchirpDecoder = CHchirpDecoder(loraSet);
% 读取文件夹下所有采样值文件
fileDir = 'E:\Pyramid_samples\Samples_base\BW125000\sf10\';
fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
writeFileDir = 'E:\Pyramid_samples\clipSignal\SF10_BW125\';
for fileCount = 1:length(fileIn)
% for fileCount = 59:59
    % 从文件中读取信号流
    [signal] = readSignalFile(fileDir, fileIn(fileCount));
    % 直接裁剪信号，因为信号已经对齐
    writeFileName = writeFileDir + "signal" + fileCount + ".sigmf-data";
    write_signal_to_file(signal(1:(12.25+20)*loraSet.dine), writeFileName);
end

toc;
fclose all;