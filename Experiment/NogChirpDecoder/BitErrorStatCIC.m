% 
clear;
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 所有采样信号的基本参数信息枚举
bw = 250;
sf = 10;
samplesRate = 2e6;
bitErrorLoc = zeros(1, sf);
valErrorStat = [];

% 读取配置和验证文件
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw * 1e3, samplesRate);
loraSet.payloadNum = 33; % payload数目

% 初始化decoder
CICDecoder = CICDecoder(loraSet);

% 读取文件夹下所有采样值文件
fileDir = 'd:\data\1_17indoor\FFT_jun\';
fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
true_bin = importdata(strcat('.\Config\bin\BW', string(bw), 'SF', string(sf), '.txt'))';

fileNum = numel(fileIn);
fileBinTrueRate = zeros(1, fileNum);
for file_i =  1: fileNum % fileNum
    [signal] = readSignalFile(fileDir, fileIn(file_i));
    off = randi([1, loraSet.dine], 1, 1);  % 窗口内随机off
    paddedSignal = [zeros(1, off) signal];
    try
        CICDecoder = CICDecoder.decode(paddedSignal);
    catch
        % disp(['文件 ' , num2str(file_i), '/', num2str(fileNum), ' 出现错误']);
        continue;
    end
    % 统计比特错误模式
    pktNum = length(CICDecoder.binRecord);
    if pktNum
        disp(['文件 ', num2str(file_i) , '/', num2str(fileNum)]);
        for binResultIndex = 1 : pktNum
            errorDist = calculateBitError(true_bin, CICDecoder.binRecord{binResultIndex}, sf);
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
% [938 14 1022 2 954 130 302 182 794 462 779 559 990 158 558 922 896 1004 839 1018 35 595 493 981 749 151 672 503 788 952 260 846 992]

toc;
fclose all;
