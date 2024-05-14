% 测试重叠信号解码
clear;
fclose all;     % 关闭所有 matlab 打开的文件
tic;            % 打开计时器

% 基本参数设置
sf = 10;
bw = 250;
samplesRate = 2e6;
bitErrorLoc = zeros(1, sf);
valErrorStat = [];
debugPath = "terminal";
DebugLevel = 4;
DebugUtil = DebugUtil(DebugLevel, debugPath);

% 读取配置和验证文件`
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw * 1e3, samplesRate);
loraSet.payloadNum = 33;  % payload数目

% 初始化decoder
obj = VarCutChirpDecoder(loraSet, DebugUtil);

% 读取文件夹下所有采样值文件
fileDir = 'd:\data\1_17indoor\FFT_jun\';
fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
true_bin = importdata(strcat('.\Config\bin\BW', string(bw), 'SF', string(sf), '.txt'))';

% 从文件中读取信号流
fileNum = numel(fileIn);
fileBinTrueRate = zeros(1, fileNum);
for file_i =  79: 79 % fileNum
    [signal] = readSignalFile(fileDir, fileIn(file_i));
    try
        obj = obj.decode(signal);
    catch
        continue;
    end
    % 统计比特错误模式
    pktNum = length(obj.payloadBin);
    if pktNum
        disp(['文件 ', num2str(file_i) , '/', num2str(fileNum)]);
        for binResultIndex = 1 : pktNum
            errorDist = calculateBitError(true_bin, cell2mat(obj.payloadBin{binResultIndex}), sf);
            if ~isempty(errorDist)
                valDist = errorDist(1, :);
                binDist = errorDist(2, :);
                disp(['valDist: ', num2str(valDist)]);
                disp(['binDist: ', num2str(binDist)]);
                valErrorStat = [valErrorStat, valDist];
                binVal = dec2bin(binDist, sf) - '0';
                bitErrorLoc = bitErrorLoc + sum(binVal, 1);
            end
        end
    % else
    %     disp(['文件 ', num2str(file_i) , '/', num2str(fileNum), ' 未检测到信号']);
    end
end
disp(bitErrorLoc);
save('errorStatistics.mat', 'valErrorStat', 'bitErrorLoc');
toc;
fclose all;
