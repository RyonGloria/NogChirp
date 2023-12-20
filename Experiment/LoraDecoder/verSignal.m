% 测试CHchirpDecoder
clear;
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 读取配置和验证文件
sf = 9;
bw = 250e3;
samplesRate = 2e6;
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
loraSet.channelNum = 2; % 信道数目
loraSet.subchirpNum = 2; % subchirp数目
loraSet.payloadNum = 18; % payload数目
LoraDecoderTmp = LoraDecoderTmp(loraSet);
true_bin = importdata(strcat('.\Code\Config\Pyramid\SF', string(sf), '.txt'))';
% 读取文件夹下所有采样值文件
fileDir = "D:\CHchirp_20Nodes\Pyramid_SF9_BW250K\Node";

for NodeIndex = 1:20
    filePath = fileDir + string(NodeIndex) + "\";
    fileIn = dir(fullfile(filePath, '*.sigmf-data'));
    for fileCount = 1:length(fileIn)
        % 从文件中读取信号流
        [signal] = readSignalFile(filePath, fileIn(fileCount));
        % 写入读取区
        LoraDecoderTmp = LoraDecoderTmp.decode(signal);

%         disp(LoraDecoderTmp.binRecord');
    
        trueRate = calculateAccuracy(true_bin, LoraDecoderTmp.binRecord);
    
        if trueRate == 1
%             disp("成功: " + filePath + fileIn(fileCount).name + " 成功率: " + trueRate);
        else
            disp("失败: " + filePath + fileIn(fileCount).name + " 成功率: " + trueRate);
        end
        
    end
end

toc;
fclose all;