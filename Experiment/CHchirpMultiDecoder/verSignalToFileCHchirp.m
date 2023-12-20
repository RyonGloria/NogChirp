% 测试CHchirpDecoder
clear;
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 配置集
% argsArr = [9, 125e3, 2, 2; ...
argsArr = [9, 125e3, 2, 4; ...
           9, 125e3, 2, 6; ...
%            9, 250e3, 2, 2; ...
           9, 250e3, 2, 4; ...
%            10, 125e3, 2, 2; ...
           10, 125e3, 2, 4; ...
           10, 125e3, 2, 6; ...
           10, 125e3, 4, 2; ...
           10, 125e3, 4, 4; ...
           10, 125e3, 4, 6; ...
           10, 250e3, 2, 4; ...
           10, 250e3, 4, 4;];
originDirPath = "D:\CHchirp_20Nodes\";
verDirPath = "D:\CHchirp_20Nodes\verSignal\";

% 遍历配置集中的所有配置
for argsIndex = 1:length(argsArr)
    sf = argsArr(argsIndex, 1);
    bw = argsArr(argsIndex, 2);
    samplesRate = 2e6;
    [loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
    loraSet.channelNum = argsArr(argsIndex, 4); % 信道数目
    loraSet.subchirpNum = argsArr(argsIndex, 3); % subchirp数目
    loraSet.payloadNum = 20; % payload数目
    DEBUG = false;
    % 从文件中加载bin groundtrurh
    true_bin = importdata(strcat('.\Code\Config\bin\SF', string(sf), '.txt'))';
    % 获取存放验证信号的目录
    if ~exist(verDirPath, 'dir')
        mkdir(verDirPath);
    end
    % 读取文件夹下对应配置的目录
    argsPath = "SF" + string(loraSet.sf) ...
                + '_BW' + string(loraSet.bw/1000) ...
                + 'K_sub' + string(loraSet.subchirpNum) ...
                + '_channels' + string(loraSet.channelNum) ...
                + '\';
    originArgsDir = originDirPath + argsPath;
    verArgsDir = verDirPath + argsPath;
    if ~exist(verArgsDir, 'dir')
        mkdir(verArgsDir);
    end
    % 初始化解调解码器
    clear CHchirpMultiDecoder;
    CHchirpMultiDecoder = CHchirpMultiDecoder(loraSet, DEBUG);
    % 遍历原始文件夹下面的所有Node节点目录
    originFileIn = dir(originArgsDir);
    for nodeIndex = 3:length(originFileIn)
        nodeName = originFileIn(nodeIndex).name;
        filePath = originArgsDir + nodeName + "\";
        % 查看验证路径是否存在
        verPath = verArgsDir + nodeName + "\";
        if ~exist(verPath, 'dir')
            mkdir(verPath);
        end
        fileIn = dir(fullfile(filePath, "*.sigmf-data"));
        for fileIndex = 1:length(fileIn) 
            % 从文件中读取信号流
            [signal] = readSignalFile(filePath, fileIn(fileIndex));
            CHchirpMultiDecoder = CHchirpMultiDecoder.decode(signal);
            % 验证解调出来的所有bin是否正确
            trueRateMax = 0;
            for binResultIndex = 1:length(CHchirpMultiDecoder.binArray)
                [true_chirp, true_rate] = vertify_bin(CHchirpMultiDecoder.binArray{binResultIndex}, true_bin);
                trueRateMax = max(trueRateMax, true_rate);
            end
            disp("文件路径:" + filePath + fileIn(fileIndex).name  + ", 正确率为：" + trueRateMax*100 + "%");
            % 验证成功
            if(trueRateMax == 1)
                fileName = verPath + string(fileIndex) + ".sigmf-data";
                write_signal_to_file(signal, fileName);
            end
        end
    end
end


toc;
fclose all;