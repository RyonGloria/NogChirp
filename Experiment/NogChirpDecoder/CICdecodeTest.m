% 
clear;
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 所有采样信号的基本参数信息枚举
bw = 250e3;
sf = 10;
samplesRate = 2e6;
% 读取配置和验证文件
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
loraSet.payloadNum = 33; % payload数目

% 初始化decoder
CICDecoder = CICDecoder(loraSet);

% fileDir = 'd:\data\ChNum_2_m2h3\';
% fileDir = 'd:\data\Collision-2_CH-2\';

fileDir = 'd:\data\1_17indoor\FFT_jun\';
fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
[signal] = readSignalFile(fileDir, fileIn(3));

off = randi([1, loraSet.dine], 1, 1);  % 窗口内随机off
disp(off);
paddedSignal = [zeros(1, off) signal];

% [938 14 1022 2 954 130 302 182 794 462 779 559 990 158 561 922 896 1004 839 571 35 595 493 981 749 151 672 503 788 952 260 846 992 ]

CICDecoder = CICDecoder.decode(paddedSignal);

dimensions = size(CICDecoder.binRecord);
lengths = cellfun(@length, CICDecoder.binRecord);
disp("⭐payloadBin dimensions: " + num2str(dimensions(1))+ "x" + num2str(dimensions(2)));
disp("⭐payloadBin lengths: " + num2str(lengths));

for i = 1:numel(CICDecoder.binRecord)
    fprintf('\n⭐Bin Cell -%d-\n', i);
    % 使用 cellfun 将当前行的每个元素格式化为字符串，并连接起来
    row_str = cellfun(@(x) sprintf('%s', mat2str(x)), CICDecoder.binRecord(i), 'UniformOutput', false);
    % 使用 strjoin 将格式化后的字符串连接起来，并输出
    disp(strjoin(row_str, ' '));
    row_str = replace(row_str, '[', '');
    row_str = replace(row_str, ']', '');
end

toc;
fclose all;
