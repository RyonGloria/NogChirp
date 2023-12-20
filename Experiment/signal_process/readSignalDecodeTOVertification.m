% 将输入文件夹中的信号文件读取分信道提取出来写入目的文件夹中
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 采样值文件读取路径和保存路径
inDir = 'D:\CHchirp_IPSN_samples\';
inDirVer = 'D:\CHchirp_IPSN_samples\';
% 设置preamble信道可选项
preambleChannelChoice = [-281.250e3, -93.750e3, 93.750e3, 281.250e3];
% 设置实验参数
PreambleChannel = 4;
payloadStartOffset = 13.25;
subchirpNum = 16;  % subchirp数目
payloadNum = 40;
channelNum = 4;
channelChoiceIndexTable = [1,2];
channelChoiceNum = length(channelChoiceIndexTable);  % 实际选择的信道数目
deleteFlag = true;
% 根据参数修改文件路径
inDir = strcat(inDir, 'sf', string(loraSet.sf), '\channel', string(channelChoiceNum), '\preamble', string(PreambleChannel), '\subchirpNum', string(subchirpNum), '\');
inDirVer = strcat(inDirVer, 'sf', string(loraSet.sf), '\channel', string(channelChoiceNum), '\preamble', string(PreambleChannel), '\subchirpNum', string(subchirpNum), '\');
% 读取配置和验证文件
[loraSet] = readLoraSet('sf7_BW125.json');

% 生成0频段的idealchirp
[downchirp, upchirp] = buildIdealchirp(loraSet, 0);
% 读取文件夹下所有采样值,bin,channel文件
fileIn = dir(fullfile(inDir, '*.sigmf-data'));     
binTxt = dir(fullfile(inDirVer, 'bin_*.txt')); 
channelTxt = dir(fullfile(inDirVer, 'channel_*.txt'));
cfoTxt = dir(fullfile(inDirVer, 'cfo_*.txt'));
% 对所有文件按创建时间排序
[resultSort] = sortFileByTime([fileIn, binTxt, channelTxt, cfoTxt]);

for fileCount = 1:length(fileIn)
    if mod(fileCount, 100) == 0           % 每循环100次，输出当前进度
        fprintf("The time is %d\n", fileCount);
    end
    % 从文件中读取信号
    [signal] = readSignalFile(inDir, fileIn(resultSort(1, fileCount)));
    % 将信号按照信道存放规则划分
    signal = reshape(signal, [], 4).';
    % 调整信号的cfo和timeoffset
    cfo = load(strcat(inDirVer, cfoTxt(resultSort(4, fileCount)).name));

    % 根据cfo生成对应0频段的idealchirp
    [downchirpCfo, upchirpCfo] = rebuild_idealchirp_cfo(loraSet, cfo, 0);
    % 获得downchirp bin
    [downchirpSync] = getDownchirpSync(loraSet, signal(PreambleChannel,:), upchirpCfo);
    % 读取对应bin.txt下的验证向量
    verBin = load(strcat(inDirVer, binTxt(resultSort(2, fileCount)).name));
    fprintf('The pkg %d bin is: ',fileCount);
    fprintf(' %d ', downchirpSync - verBin(1));
    % 获得downchirpbin对应的信道向量
    [channelChoice] = createChannelChoiceVector(downchirpSync);
    % 解bin
    [chirp1Bin] = decodeSubchirp(loraSet, signal, channelChoice, channelChoiceIndexTable, downchirpCfo, payloadStartOffset, payloadNum, subchirpNum, channelChoiceNum);
    fprintf(' %d ', chirp1Bin - verBin);
    binTrueNum = sum((chirp1Bin - verBin) == 0);
    % 识别channel
    verChannel = load(strcat(inDirVer, channelTxt(resultSort(3, fileCount)).name));
    fprintf('channel is: ');
    [channelPredictArray] = predictChannelVertex(loraSet, signal, downchirpCfo, payloadStartOffset, payloadNum, subchirpNum);
    for channelCount = 1:size(channelPredictArray, 1)
        fprintf(' %d ', channelPredictArray(channelCount, :) - verChannel(channelCount, :));
    end
    channelTrueNum = sum((channelPredictArray - verChannel) == 0, "all");
    fprintf('\n');

    if (binTrueNum ~= payloadNum || channelTrueNum ~= subchirpNum*payloadNum) && deleteFlag == true
        delete(strcat(inDir, fileIn(resultSort(1, fileCount)).name));
        delete(strcat(inDirVer, binTxt(resultSort(2, fileCount)).name));
        delete(strcat(inDirVer, channelTxt(resultSort(3, fileCount)).name));
        delete(strcat(inDirVer, cfoTxt(resultSort(4, fileCount)).name));
    end

end
fclose all;