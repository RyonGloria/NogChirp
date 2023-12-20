% 测试CHchirpMultiDecoder解单包的性能，保存实验结果
clear;
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 参数枚举
argsArr = [
%     9, 125e3, 2, 2; ...
%     9, 125e3, 2, 4; ...
%     9, 125e3, 2, 6; ...
%     9, 250e3, 2, 2; ...
%     9, 250e3, 2, 4; ...
    10, 125e3, 2, 2; ...
    % 10, 125e3, 2, 4; ...
%     10, 125e3, 2, 6; ...
%     10, 125e3, 4, 2; ...
    % 10, 125e3, 4, 4; ...
%     10, 125e3, 4, 6; ...
%     10, 250e3, 2, 4; ...
%     10, 250e3, 4, 4;
    ];

dirPath = "D:\CHchirp_20Nodes\verSignal\SF";
savePath = ".\Result\detectExp\";
debugPath = "terminal";
round = 40;
collNum = 8; % 冲突数目固定为8
DebugLevel = 4;
DebugUtil = DebugUtil(DebugLevel, debugPath);
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

% 枚举参数实验
for argsIndex = 1:size(argsArr, 1)
    % 所有采样信号的基本参数信息枚举
    sf = argsArr(argsIndex, 1);
    bw = argsArr(argsIndex, 2);
    subchirpNum = argsArr(argsIndex, 3);
    channelNum = argsArr(argsIndex, 4);
    samplesRate = 2e6;
    
    
    record = cell(1, 0);  % 记录每一次实验的过程值
    SNR = -10;
    result = cell(1, 0);  % 记录每轮的结果
    
    % 读取配置和验证文件
    [loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
    loraSet.channelNum = channelNum; % 信道数目
    loraSet.subchirpNum = subchirpNum; % subchirp数目
    loraSet.payloadNum = 20; % payload数目
    clear CHchirpMultiDecoder;
    CHchirpMultiDecoder = CHchirpMultiDecoder(loraSet, DebugUtil);
    SignalLength = loraSet.dine*80;  % 整个信号的最大长度
    argsName = 'SF' + string(loraSet.sf) ...
        + '_BW' + string(loraSet.bw/1000) ...
        + '_subchirp' + string(loraSet.subchirpNum) ...
        + '_channel' + string(loraSet.channelNum) ...
        + '_CollNum' + string(collNum);
    % 从文件中加载bin groundtrurh
    true_bin = importdata(strcat('.\Code\Config\bin\SF', string(sf), '.txt'))';
    % 读取文件夹下所有采样值文件
    fileDir = dirPath + string(loraSet.sf) ...
        + '_BW' + string(loraSet.bw/1000) ...
        + 'K_sub' + string(loraSet.subchirpNum) ...
        + '_channels' + string(loraSet.channelNum);
    if exist(fileDir, "dir")
        fileDir = fileDir + "\";
        fileIn = dir(fullfile(fileDir));
        count = 0;  % 记录每个冲突情况的次数
        countAll = 0; % 记录这次实验的所有次数
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
            for time = 1:round
                count = count + 1;
                countAll = countAll + 1;
                % 合成冲突
                NodeIndex = randi([3, length(fileIn)], 1, collNum);
                fileIndex = zeros(1, collNum);
                % fileIndex = randi([1, length(fileIn)], 1, collNum);
                winoff = randi([1, 30], 1, collNum); % 随机窗口
                off = randi([1, loraSet.dine], 1, collNum); % 窗口内随机off
                AMP = zeros(1, collNum);
                SIR = (rand(1, collNum) * 10) - 5;
                signalAll = zeros(1, SignalLength);
                % 从文件中随机读取pkgNum个信号,设置随机偏移，合成冲突信号
                for pkgIndex = 1:collNum
                    fileDirPath = fileDir + fileIn(NodeIndex(pkgIndex)).name + "\";
                    fileFull = dir(fullfile(fileDirPath, '*.sigmf-data'));
                    fileRandiIndex = randi([1, length(fileFull)], 1, 1);
                    fileIndex(1, pkgIndex) = fileRandiIndex;
                    signal = readSignalFile(fileDirPath, fileFull(fileIndex(pkgIndex)));
                    signal = signal(1:loraSet.dine*40);
                    if pkgIndex == 1
                        % 添加噪声
                        mainSignalAmp = mean(abs(signal(3*loraSet.dine : 5*loraSet.dine)));
                        noiseAmp = mainSignalAmp/10^(SNR/20);
                        noise = (noiseAmp/sqrt(2) * randn([1 SignalLength]) + 1i*noiseAmp/sqrt(2) * randn([1 SignalLength]));   % 生成噪声
                        offset = loraSet.dine*winoff(pkgIndex) + off(pkgIndex);
                        signal = [zeros(1, offset), signal, zeros(1, SignalLength - length(signal) - offset)];
                        signal = signal + noise;
                    else % 不是第一个包，需要根据SIR设置强度
                        signalAmp = mean(abs(signal(3*loraSet.dine : 5*loraSet.dine)));
                        amp = mainSignalAmp/(10^(SIR(pkgIndex)/20))/signalAmp;
                        AMP(pkgIndex) = amp;
                        signal = amp*signal;
                        offset = loraSet.dine*winoff(pkgIndex) + off(pkgIndex);
                        signal = [zeros(1, offset), signal, zeros(1, SignalLength - length(signal) - offset)];
                    end
                    signalAll = signalAll + signal;
                end
                % 解码
                try
                    % CHchirpMultiDecoder = CHchirpMultiDecoder.decode(signalAll);
                    CHchirpMultiDecoder = CHchirpMultiDecoder.decodeCalculateDetect(signalAll);
                catch
                    disp("参数：" + argsName ...
                        + "\n文件:" + fileIndex ...
                        + "\nwinoff: " + winoff ...
                        + "\noff" + off ...
                        + "\n power" + AMP + "出现错误");
                    decodeFlag = true;
                end
                
                if decodeFlag == false
                    % recordRateTmp = zeros(1, length(CHchirpMultiDecoder.binArray));
                    % for binResultIndex = 1:length(CHchirpMultiDecoder.binArray)
                    %     [true_chirp, recordRateTmp(binResultIndex)] = vertify_bin(CHchirpMultiDecoder.binArray{binResultIndex}, true_bin);
                    % end
                    % recordRateTmp = sort(recordRateTmp, 'descend');
                    % trueRate = 0;
                    % for pkgIndex = 1:collNum
                    %     if length(recordRateTmp) >= pkgIndex
                    %         trueRate = trueRate + recordRateTmp(pkgIndex);
                    %     end
                    % end
                    % trueRate = trueRate / collNum;
                    if CHchirpMultiDecoder.detectSFDCount >= collNum
                        trueRate = 1;
                    else
                        trueRate = CHchirpMultiDecoder.detectSFDCount/collNum;
                    end
                    % detectNum = CHchirpMultiDecoder.detectCount;
                else
                    % recordRateTmp = 0;
                    trueRate = 0;
                    % detectNum = 0;
                    countError = countError + 1;
                    countErrorAll = countErrorAll + 1;
                    countErrorTmp = countErrorTmp + 1;
                end
                
                % 记录结果
                record{1, countAll} = SNR;
                record{2, countAll} = argsName;
                record{3, countAll} = NodeIndex;
                record{4, countAll} = fileIndex;
                record{5, countAll} = winoff;
                record{6, countAll} = off;
                record{7, countAll} = collNum;
                record{8, countAll} = SIR;
                record{9, countAll} = AMP;
                % record{10, countAll} = detectNum;
                % record{10, countAll} = CHchirpMultiDecoder.binArray;
                record{10, countAll} = trueRate;
                % record{12, countAll} = recordRateTmp;
                record{11, countAll} = decodeFlag;
                decodeFlag = false;
                
                trueRateNow = trueRateNow + trueRate;
                disp(   "参数: " + argsName ...
                    + ", SNR: " + SNR ...
                    + ", 总轮次计数: " + countAll ...
                    + ", 当前轮次计数: " + count ...
                    + ", 所有轮次错误次数：" + countErrorAll ...
                    + ", 当前轮次错误次数：" + countError ...
                    + ", 冲突发现正确率: " + trueRate);
            end
            % 判断是否进行下一轮，或者直接退出程序
            
            % 趋于稳定
            if count >= 40 && abs(trueRateAll/((count-countError)-(round-countErrorTmp)) - (trueRateAll+trueRateNow)/(count-countError)) < 0.1
                convergeCount = convergeCount + 1;
                if convergeCount >= 3
                    result{SNR + 10 + 1} = (trueRateAll+trueRateNow)/(count-countError);
                    save(savePath + datestr(now, 30) + argsName + ".mat", 'record', "result");
                    disp("当前冲突数目结果收敛: " + result{SNR + 10 + 1});
                    % SNR超过10，实验结束
                    if SNR >= 10
                        disp("当前冲突数目大于等于14,实验结束");
                        break;
                    end
                    count = 0;
                    countError = 0;
                    trueRateAll = 0;
                    SNR = SNR + 5;
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
            plotResult(i) = result{i};
        end
        plot(plotResult);
        title(argsName+'成功率分布');
        xlabel('冲突数目');
        ylabel('symbol正确率');
        saveas(gcf, savePath + datestr(now, 30) + argsName + ".fig");
        saveas(gcf, savePath + datestr(now, 30) + argsName + ".png");
        save(savePath + datestr(now, 30) + argsName + ".mat", "record", "result");
    end
end

toc;
fclose all;