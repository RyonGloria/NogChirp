% 将输入文件夹中的信号文件读取分信道提取出来写入目的文件夹中
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 采样值文件读取路径和保存路径
inDir = 'E:\share\';
inDirVer = 'E:\uhd_share\';
writeDir = 'D:\CHchirp_IPSN_samples\';
% 设置preamble信道可选项
% preambleChannelChoice = [175e3, 375e3, 575e3, 775e3];
preambleChannelChoice = [-300e3, -100e3, 100e3, 300e3];
% 设置实验参数
preambleChannel = 2;
payloadStartOffset = 13.25;
subchirpNum = 1;
payloadNum = 40;
channelChoiceNum = 4;
% channelChoiceIndexTable = [1,2,3,4];
channelChoiceIndexTable = [2,2,2,2];
writeFlag = false;  % 是否进行写入标志位
% 读取配置和验证文件
[loraSet] = readLoraSet('sf7_BW125.json');

% 生成0频段的idealchirp
[downchirp, upchirp] = buildIdealchirp(loraSet, 0);
% 读取文件夹下所有采样值,bin,channel文件
fileIn = dir(fullfile(inDir, '*.sigmf-data'));     
binTxt = dir(fullfile(inDirVer, 'bin_*.txt')); 
channelTxt = dir(fullfile(inDirVer, 'channel_*.txt')); 
% 对所有文件按创建时间排序
[resultSort] = sortFileByTime([fileIn, binTxt, channelTxt]);

for fileCount = 1:1
    if mod(fileCount, 100) == 0           % 每循环100次，输出当前进度
        fprintf("The time is %d\n", fileCount);
    end
    % 从文件中读取信号文件
    [signal] = readSignalFile(inDir, fileIn(resultSort(1, fileCount)));
    % 调整信号的cfo和timeoffset
    [cfo, signal, packageFlag] = alignSignal(loraSet, signal, downchirp, upchirp, preambleChannelChoice(preambleChannel));
%     if packageFlag == false
%         continue;
%     end
    stft(signal(1:40*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
    fs = 2e6;
    N = loraSet.dine;
    n = 0:N-1;
    t = 0:1/fs:1-1/fs;
    f = n*fs/N;
    b = fir1(30, 0.05, "low");
    singalOut = filter(b, 1, signal);
    stft(singalOut(1:40*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
end
fclose all;