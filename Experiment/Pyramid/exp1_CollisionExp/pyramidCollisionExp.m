% 
clear;
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 所有采样信号的基本参数信息枚举
bw = 125e3;
sf = 10;
channelNum = 2;
samplesRate = 2e6;
DEBUG = false;
record = cell(3, 1);  % 设置元胞数组记录结果
SNR = 5;
result = cell(1, 0);

% 读取配置和验证文件
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
loraSet.channelNum = channelNum; % 信道数目
loraSet.payloadNum = 20; % payload数目
SignalLength = loraSet.dine*80;  % 整个信号的最大长度
argsName = 'SF' + string(loraSet.sf) ...
    + '_BW' + string(loraSet.bw/1000) ...
    + 'CHN' + string(channelNum);

PyramidDecoder = PyramidDecoder(loraSet);
% 从文件中加载bin groundtrurh
true_bin = importdata(strcat('.\Code\Config\bin\SF', string(sf), '.txt'))';
% 读取文件夹下所有采样值文件
fileDir = 'D:\VMwareCopy\verSignal\SF' + string(loraSet.sf) ...
    + 'BW' + string(loraSet.bw/1000);
if exist(fileDir, "dir")
    fileDir = fileDir + "\";
    fileIn = dir(fullfile(fileDir));
    count = 0;  % 记录每个冲突情况的次数
    countAll = 0; % 记录这次实验的所有次数
    collNum = 1; % 记录当前轮次冲突的数目
    countError = 0;
    countErrorAll = 0;
    trueRateAll = 0;
    convergeCount = 0;
    decodeFlag = false;

    % 开始循环进行冲突的测量
    while true
        trueRateNow = 0;
        countErrorTmp = 0;
        % 每轮运行10次，看结果是否继续
        for time = 1:10
            count = count + 1;
            countAll = countAll + 1;
            % 设置随机信道
            preambleChannel = randi([1, channelNum], 1, collNum); % 随机设置信号所在信道的位置
            preambleChannelNum = zeros(1, channelNum);
            % 统计每个信道的信号数目
            for channelIndex = 1:collNum
                preambleChannelNum(preambleChannel(channelIndex)) = preambleChannelNum(preambleChannel(channelIndex)) + 1;
            end
            % 对每个信道进行单独解码，统计结果
            channelBinTrueRate = zeros(1, channelNum);
            % 记录参数
            fileIndexRecord = cell(1, channelNum);
            winoffRecord = cell(1, channelNum);
            offRecord = cell(1, channelNum);
            SIRRecord = cell(1, channelNum);
            SNRrecord = cell(1, channelNum);
            binRecord = cell(1, channelNum);
            SRRRecord = cell(1, channelNum);
            for channelIndex = 1:channelNum
                if preambleChannelNum(channelIndex) > 0
                    % 设置随机数
                    NodeIndex = randi([3, length(fileIn)], 1, preambleChannelNum(channelIndex));
                    fileIndex = zeros(1, preambleChannelNum(channelIndex));
                    winoff = randi([2, 30], 1, preambleChannelNum(channelIndex)); % 随机窗口
                    off = randi([1, loraSet.dine], 1, preambleChannelNum(channelIndex)); % 窗口内随机off
                    SIR = (rand(1, preambleChannelNum(channelIndex)) * 10) - 5;
                    % 合成信号
                    signalAll = zeros(1, SignalLength);
                    % 从文件中随机读取pkgNum个信号,设置随机偏移，合成冲突信号
                    for pkgIndex = 1:preambleChannelNum(channelIndex)
                        fileDirPath = fileDir + fileIn(NodeIndex(pkgIndex)).name + "\";
                        fileFull = dir(fullfile(fileDirPath, '*.sigmf-data'));
                        fileRandiIndex = randi([1, length(fileFull)], 1, 1);
                        fileIndex(1, pkgIndex) = fileRandiIndex;
                        signal = readSignalFile(fileDirPath, fileFull(fileRandiIndex));
                        signal = signal(1 : loraSet.dine*40);
                        if pkgIndex == 1
                            % 添加噪声
                            mainSignalAmp = mean(abs(signal(3*loraSet.dine : 5*loraSet.dine)));
                            noiseAmp = mainSignalAmp/10^(SNR/20);
                            noise = (noiseAmp/sqrt(2) * randn([1 SignalLength]) + 1i*noiseAmp/sqrt(2) * randn([1 SignalLength]));   % 生成噪声
%                             signal = signal + noise;
                            offset = loraSet.dine*winoff(pkgIndex) + off(pkgIndex);
                            signal = [zeros(1, offset), signal, zeros(1, SignalLength - length(signal) - offset)];
                            signal = signal + noise;
                        else % 不是第一个包，需要根据SIR设置强度
                            signalAmp = mean(abs(signal(3*loraSet.dine : 5*loraSet.dine)));
                            amp = mainSignalAmp/(10^(SIR(pkgIndex)/20))/signalAmp;
                            signal = amp*signal;
                            offset = loraSet.dine*winoff(pkgIndex) + off(pkgIndex);
                            signal = [zeros(1, offset), signal, zeros(1, SignalLength - length(signal) - offset)];
                        end
                        signalAll = signalAll + signal;
                    end
                    % 解码
                    try
                        PyramidDecoder = PyramidDecoder.decode(signalAll);
                    catch
                        disp("参数：" + argsName ...
                            + "\n文件:" + fileIndex ...
                            + "\nwinoff: " + winoff ...
                            + "\noff" + off ...
                            + "\nchannel" + preambleChannel ...
                            + "出现错误");
                        decodeFlag = true;
                    end
                    % 计算该信道的正确率
                    recordRateTmp = zeros(1, length(PyramidDecoder.binRecord));
                    for binResultIndex = 1:length(PyramidDecoder.binRecord)
                        recordRateTmp(binResultIndex) = calculateAccuracy(true_bin, PyramidDecoder.binRecord{binResultIndex});
                    end
                    recordRateTmp = sort(recordRateTmp, 'descend');
                    trueRate = 0;
                    for pkgIndex = 1:preambleChannelNum(channelIndex)
                        if length(recordRateTmp) >= pkgIndex
                            trueRate = trueRate + recordRateTmp(pkgIndex);
                        end
                    end
                    trueRate = trueRate / preambleChannelNum(channelIndex);
                    channelBinTrueRate(channelIndex) = trueRate;
                    % 记录当前信道信息
                    fileIndexRecord{channelIndex} = fileIndex;
                    winoffRecord{channelIndex} = winoff;
                    offRecord{channelIndex} = off;
                    SIRRecord{channelIndex} = SIR;
                    SNRrecord{channelIndex} = SNR;
                    binRecord{channelIndex} = PyramidDecoder.binRecord;
                    SRRRecord{channelIndex} = recordRateTmp;
                else
                    % 记录当前信道信息
                    trueRate = 0;
                    fileIndexRecord{channelIndex} = 0;
                    winoffRecord{channelIndex} = 0;
                    offRecord{channelIndex} = 0;
                    SIRRecord{channelIndex} = 0;
                    SNRrecord{channelIndex} = 0;
                    binRecord{channelIndex} = 0;
                    channelBinTrueRate(channelIndex) = 0;
                    SRRRecord{channelIndex} = 0;
                end
            end
            % 计算总正确率
            trueRateRoundAll = 0;
            if decodeFlag == false
                for channelIndex = 1:channelNum
                    trueRateRoundAll = trueRateRoundAll + channelBinTrueRate(channelIndex)*(preambleChannelNum(channelIndex)/collNum);
                end
            else
                trueRateRoundAll = 0;
                countError = countError + 1;
                countErrorAll = countErrorAll + 1;
                countErrorTmp = countErrorTmp + 1;
                decodeFlag = false;
            end
            
            % 记录结果
            record{1, countAll} = collNum;
            record{2, countAll} = argsName;
            record{3, countAll} = fileIndexRecord;
            record{4, countAll} = winoffRecord;
            record{5, countAll} = offRecord;
            record{7, countAll} = SIRRecord;
            record{6, countAll} = SNRrecord;
            record{8, countAll} = binRecord;
            record{9, countAll} = trueRate;
            record{10, countAll} = SRRRecord;
            record{11, countAll} = decodeFlag;

            trueRateNow = trueRateNow + trueRateRoundAll;
            disp(   "参数: " + argsName ...
                    + ", 冲突数目: " + collNum ...
                    + ", 总轮次计数: " + countAll ...
                    + ", 当前轮次计数: " + count ...
                    + ", 所有轮次错误次数：" + countErrorAll ...
                    + ", 当前轮次错误次数：" + countError ...
                    + ", 成功率: " + trueRateRoundAll);
        end
        % 判断是否进行下一轮，或者直接退出程序

        % 趋于稳定
        if count >= 40 && abs(trueRateAll/((count-countError)-(10-countErrorTmp)) - (trueRateAll+trueRateNow)/(count-countError)) < 0.01
            convergeCount = convergeCount + 1;
            if convergeCount >= 3
                result{collNum} = (trueRateAll+trueRateNow)/(count-countError);
                save(".\Result\PyramidCollisionExp\" + datestr(now, 30) + argsName + ".mat", 'record', "result");
                disp("当前冲突数目结果收敛: " + result);
                % 正确率低于阈值，结束程序
                if (trueRateAll+trueRateNow)/(count-countError) < 0.4
                    disp("当前正确率低于0.4,实验结束");
                    break;
                end
                count = 0;
                countError = 0;
                trueRateAll = 0;
                collNum = collNum + 1;
                convergeCount = 0;
            else
                trueRateAll = trueRateAll + trueRateNow;
            end
        else
            trueRateAll = trueRateAll + trueRateNow;
            convergeCount = 0;
        end
        
    end
    plotResult = zeros(1, length(result));
    for i = 1:length(result)
        plotResult = result{i};
    end
    plot(plotResult);
    title(argsName+'成功率分布');
    xlabel('冲突数目');
    ylabel('symbol正确率');
    saveas(gcf, ".\Result\PyramidCollisionExp\" + datestr(now, 30) + argsName + ".fig");
    saveas(gcf, ".\Result\PyramidCollisionExp\" + datestr(now, 30) + argsName + ".png");
    save(".\Result\PyramidCollisionExp\" + datestr(now, 30) + argsName + ".mat", "record", "result");
end

toc;
fclose all;