% 将输入文件夹中的信号文件读取分信道提取出来写入目的文件夹中
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 采样值文件读取路径和保存路径
inDir = 'E:\share\';
inDirVer = 'E:\uhd_share\';
writeDir = 'D:\CHchirp_IPSN_samples\';
% 设置preamble信道可选项
preambleChannelChoice = [-281.250e3, -93.750e3, 93.750e3, 281.250e3];
% 设置实验参数
preambleChannel = 4;
payloadStartOffset = 13.25;
subchirpNum = 16;
payloadNum = 40;
channelNum = 4;
channelChoiceIndexTable = [1,2];
channelChoiceNum = length(channelChoiceIndexTable);
writeFlag = true;  % 是否进行写入标志位
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

for fileCount = 1:length(fileIn)
    if mod(fileCount, 100) == 0           % 每循环100次，输出当前进度
        fprintf("The time is %d\n", fileCount);
    end
    % 从文件中读取信号文件
    [signal] = readSignalFile(inDir, fileIn(resultSort(1, fileCount)));
    % 调整信号的cfo和timeoffset
    [cfo, signal, packageFlag] = alignSignal(loraSet, signal, downchirp, upchirp, preambleChannelChoice(channelNum-preambleChannel+1));
    if packageFlag == false
        continue;
    end
    % 将信号根据信道滤波划分
    [singalOut] = divideChannel(loraSet, signal, preambleChannelChoice, false);
    % 根据cfo生成对应0频段的idealchirp
    [downchirpCfo, upchirpCfo] = rebuild_idealchirp_cfo(loraSet, cfo, 0);
    % 获得downchirp bin
    [downchirpSync] = getDownchirpSync(loraSet, singalOut(preambleChannel,:), upchirpCfo);
    % 读取对应bin.txt下的验证向量
    verBin = load(strcat(inDirVer, binTxt(resultSort(2, fileCount)).name));
    fprintf('The pkg %d bin is: ',fileCount);
    fprintf(' %d ', downchirpSync - verBin(1));
    % 获得downchirpbin对应的信道向量
    [channelChoice] = createChannelChoiceVector(downchirpSync);
    % 解bin
    [chirp1Bin] = decodeSubchirp(loraSet, singalOut, channelChoice, channelChoiceIndexTable, downchirpCfo, payloadStartOffset, payloadNum, subchirpNum, channelChoiceNum);
    fprintf(' %d ', chirp1Bin - verBin');
    binTrueNum = sum((chirp1Bin - verBin') == 0);
    % 识别channel
    verChannel = load(strcat(inDirVer, channelTxt(resultSort(3, fileCount)).name));
    verChannel(:) = channelChoiceIndexTable(verChannel(:)+1);
    fprintf('channel is: ');
    [channelPredictArray] = predictChannelVertex(loraSet, singalOut, downchirpCfo, payloadStartOffset, payloadNum, subchirpNum);
    for channelCount = 1:size(channelPredictArray, 1)
        fprintf(' %d ', channelPredictArray(channelCount, :) - verChannel(channelCount, :));
    end
    channelTrueNum = sum((channelPredictArray - verChannel) == 0, "all");
    fprintf('\n');

    if binTrueNum == payloadNum && channelTrueNum == subchirpNum*payloadNum && writeFlag == true
        signal = reshape(singalOut.', 1, []);
        writeSignalFile(loraSet, signal, writeDir, preambleChannel, channelChoiceNum, subchirpNum, downchirpSync);
        writeBinChannel(loraSet, chirp1Bin, channelPredictArray, writeDir, preambleChannel, channelChoiceNum, subchirpNum, downchirpSync, cfo);
    end

end
fclose all;