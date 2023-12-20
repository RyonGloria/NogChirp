% 将输入文件夹中的信号文件读取分信道提取出来写入目的文件夹中
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 采样值文件读取路径和保存路径
inDir = 'E:\share\samples\';
inDirVer = 'E:\uhd_share\samples\';
writeDir = 'D:\CHchirp_IPSN_samples\BW250\';
% 读取配置和验证文件
[loraSet] = readLoraSet('sf8_BW250.json');
% 设置preamble信道可选项
bw = loraSet.bw;
gbw = bw/2;
preambleChannelChoice = [-(bw+gbw)/2-(bw+gbw), -(bw+gbw)/2, (bw+gbw)/2, (bw+gbw)/2+(bw+gbw)];
% preambleChannelChoice = [-(bw+gbw)/2, (bw+gbw)/2];
% 设置实验参数
preambleChannel = 1;
payloadStartOffset = 13.25;
subchirpNum = 1;
payloadNum = 40;
channelNum = 4;
% channelNum = 2;
channelChoiceIndexTable = [1,2,3,4];
% channelChoiceIndexTable = [1,2];
channelChoiceNum = length(channelChoiceIndexTable);
% channelChoiceNum = 3;
writeFlag = true;  % 是否进行写入标志位

% 生成0频段的idealchirp
[downchirp, upchirp] = buildIdealchirp(loraSet, 0);

while true
    % 获取文件名
    checkDir = strcat(writeDir, 'sf', string(loraSet.sf), '\channel', string(channelChoiceNum), '\preamble', string(preambleChannel), '\subchirpNum', string(subchirpNum), '\');
    file = dir(fullfile(checkDir, '*.sigmf-data')); 
    numArray = zeros(1, loraSet.fft_x);
    % 构建索引表
    for i = 1:2^(loraSet.sf - 7):loraSet.fft_x
        numArray(i) = i;
    end
    % 根据索引表检查缺失文件
    for i = 1:length(file)
        tmpName = file(i).name;
        tmpNameStr = tmpName(isstrprop(tmpName, 'digit'));
        tmpName = str2num(tmpNameStr);

        binFile = strcat(checkDir,'bin_downchirpsync', tmpNameStr, '.txt');
        channelFile = strcat(checkDir,'channel_downchirpsync', tmpNameStr, '.txt');
        cfoFile = strcat(checkDir,'cfo_downchirpsync', tmpNameStr, '.txt');
        if exist(binFile, 'file') && exist(channelFile, 'file') && exist(cfoFile, 'file')
            numArray(tmpName) = 0;
        end
    end
    numArray = numArray((numArray > 0));
    numArray = numArray - 1;
    fprintf(" %d,", numArray);
    if ~isempty(numArray)
        writematrix(numArray,'E:\share\notDone.txt','Delimiter',',');
        writematrix(numArray,'E:\uhd_share\notDone.txt','Delimiter',',');
        fprintf("Matlab write Done!\n");
    else
        numArray = zeros(1, 128);
        for i = 0:127
            numArray(i+1) = i*2^(loraSet.sf - 7);
        end
        writematrix(numArray,'E:\share\Done.txt','Delimiter',',');
        writematrix(numArray,'E:\uhd_share\Done.txt','Delimiter',',');
        fprintf("Matlab write Done!\n");
        subchirpNum = subchirpNum*2;
        if(subchirpNum > loraSet.factor)
            subchirpNum = 1;
            preambleChannel = preambleChannel + 1;
            if(preambleChannel > channelChoiceNum)
                preambleChannel = 1;
                channelChoiceNum = channelChoiceNum - 1;
                channelChoiceIndexTable = channelChoiceIndexTable(1:channelChoiceNum);
                if(channelChoiceNum < 2)
                    break;
                end
            end
        end
    end
    fclose all;

    fprintf("wait gr-lora...\n");
    while ~exist('E:\share\grloraDone.txt','file')
    end
    % 等待一会保证文件写成功
    pause(1);
    delete('E:\share\grloraDone.txt');
    fprintf("gr-lora Done!\n");
    % 读取文件夹下所有采样值,bin,channel文件
    fileIn = dir(fullfile(inDir, '*.sigmf-data'));     
    binTxt = dir(fullfile(inDirVer, 'bin_*.txt')); 
    channelTxt = dir(fullfile(inDirVer, 'channel_*.txt')); 
    % 解决文件数目不对的问题
    lengthMin = min([length(fileIn), length(binTxt), length(channelTxt)]);
    fileIn = fileIn(1:lengthMin);
    binTxt = binTxt(1:lengthMin);
    channelTxt = channelTxt(1:lengthMin);
    % 对所有文件按创建时间排序
    [resultSort] = sortFileByTime([fileIn, binTxt, channelTxt]);
    
    for fileCount = 1:length(fileIn)
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
    
        if binTrueNum == payloadNum && channelTrueNum == subchirpNum*payloadNum && downchirpSync - verBin(1) == 0 && writeFlag == true
            signal = reshape(singalOut.', 1, []);
            writeSignalFile(loraSet, signal, writeDir, preambleChannel, channelChoiceNum, subchirpNum, downchirpSync);
            writeBinChannel(loraSet, chirp1Bin, channelPredictArray, writeDir, preambleChannel, channelChoiceNum, subchirpNum, downchirpSync, cfo);
        end
    
    end
    pause(1);
    fclose all;
    
    % 读取文件夹下所有采样值,bin,channel文件
    deleteFlag = true;
    writeDirTmp = strcat(writeDir, 'sf', string(loraSet.sf), '\channel', string(channelChoiceNum), '\preamble', string(preambleChannel), '\subchirpNum', string(subchirpNum), '\');
    fileIn = dir(fullfile(writeDirTmp, '*.sigmf-data'));     
    binTxt = dir(fullfile(writeDirTmp, 'bin_*.txt')); 
    channelTxt = dir(fullfile(writeDirTmp, 'channel_*.txt'));
    cfoTxt = dir(fullfile(writeDirTmp, 'cfo_*.txt'));
    % 解决文件数目不对的问题
    lengthMin = min([length(fileIn), length(binTxt), length(channelTxt), length(cfoTxt)]);
    fileIn = fileIn(1:lengthMin);
    binTxt = binTxt(1:lengthMin);
    channelTxt = channelTxt(1:lengthMin);
    cfoTxt = cfoTxt(1:lengthMin);
    % 对所有文件按创建时间排序
    [resultSort] = sortFileByTime([fileIn, binTxt, channelTxt, cfoTxt]);
    
    for fileCount = 1:length(fileIn)
        % 从文件中读取信号
        [signal] = readSignalFile(writeDirTmp, fileIn(resultSort(1, fileCount)));
        % 将信号按照信道存放规则划分
        signal = reshape(signal, [], 4).';
        % 调整信号的cfo和timeoffset
        cfo = load(strcat(writeDirTmp, cfoTxt(resultSort(4, fileCount)).name));
    
        % 根据cfo生成对应0频段的idealchirp
        [downchirpCfo, upchirpCfo] = rebuild_idealchirp_cfo(loraSet, cfo, 0);
        % 获得downchirp bin
        [downchirpSync] = getDownchirpSync(loraSet, signal(preambleChannel,:), upchirpCfo);
        % 读取对应bin.txt下的验证向量
        verBin = load(strcat(writeDirTmp, binTxt(resultSort(2, fileCount)).name));
        fprintf('The pkg %d bin is: ',fileCount);
        fprintf(' %d ', downchirpSync - verBin(1));
        % 获得downchirpbin对应的信道向量
        [channelChoice] = createChannelChoiceVector(downchirpSync);
        % 解bin
        [chirp1Bin] = decodeSubchirp(loraSet, signal, channelChoice, channelChoiceIndexTable, downchirpCfo, payloadStartOffset, payloadNum, subchirpNum, channelChoiceNum);
        fprintf(' %d ', chirp1Bin - verBin);
        binTrueNum = sum((chirp1Bin - verBin) == 0);
        % 识别channel
        verChannel = load(strcat(writeDirTmp, channelTxt(resultSort(3, fileCount)).name));
        fprintf('channel is: ');
        [channelPredictArray] = predictChannelVertex(loraSet, signal, downchirpCfo, payloadStartOffset, payloadNum, subchirpNum);
        for channelCount = 1:size(channelPredictArray, 1)
            fprintf(' %d ', channelPredictArray(channelCount, :) - verChannel(channelCount, :));
        end
        channelTrueNum = sum((channelPredictArray - verChannel) == 0, "all");
        fprintf('\n');
    
        if (binTrueNum ~= payloadNum || channelTrueNum ~= subchirpNum*payloadNum || downchirpSync - verBin(1) ~= 0) && deleteFlag == true
            delete(strcat(writeDirTmp, fileIn(resultSort(1, fileCount)).name));
            delete(strcat(writeDirTmp, binTxt(resultSort(2, fileCount)).name));
            delete(strcat(writeDirTmp, channelTxt(resultSort(3, fileCount)).name));
            delete(strcat(writeDirTmp, cfoTxt(resultSort(4, fileCount)).name));
        end
    end
    
end