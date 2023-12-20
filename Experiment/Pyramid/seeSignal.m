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
fileDir = 'E:\share\samples\';
% fileDir = 'E:\Pyramid_samples\clipSignal\SF10_BW125\';
fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
for fileCount = 1:length(fileIn)
% for fileCount = 1:1
    % 从文件中读取信号流
    [signal] = readSignalFile(fileDir, fileIn(fileCount));

%     FFT_plot(signal(12.25*loraSet.dine+1:end), loraSet, CHchirpDecoder.idealDownchirp, 20);

    figure(1);
    stft(signal(1:(12.25+24)*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);

end

toc;
fclose all;