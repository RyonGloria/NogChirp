% 
clear;
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 所有采样信号的基本参数信息枚举
bw = 125e3;
sf = 10;
samplesRate = 2e6;
round = 40;
record = cell(3, 1);  % 设置元胞数组记录结果
% SNR = 5;
result = cell(1, 0);

% 读取配置和验证文件
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
loraSet.payloadNum = 23; % payload数目
SignalLength = loraSet.dine*80;  % 整个信号的最大长度

% 初始化decoder
CICDecoder = CICDecoder(loraSet);

% 读取文件夹下所有采样值文件
% fileDir = '\\192.168.3.102\e\data\ChNum_1_m2\';
% fileDir = '\\192.168.3.102\e\data\nodelay_231219\';
% fileDir = '\\192.168.3.102\e\data\ChNum_3_l1m2h3\';
% fileDir = '\\192.168.3.102\e\data\ChNum_3_l2m3h1\';
fileDir = 'd:\data\ChNum_2_m2h3\';
fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
% 从文件中读取信号流
[signal] = readSignalFile(fileDir, fileIn(1));

emptySignal = zeros(1, 10000); % create an array of zeros with the specified length
paddedSignal = [emptySignal signal]; % concatenate the empty signal with the original signal

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
end

% [810 722 386 614 850 593 126 1022 842 677 347 862 581 551 574 113 806 827 822 85 421 587 421]
% [810 1010 386 614 850 406 126 1022 367 677 347 284 27 551 1004 1006 49 836 822 85 265 587 421]
% 准确率 14/23 20/23  0/23

% ⭐Bin Cell -1-
% [810 1010 386 614 850 406 126 1022 367 677 347 284 27 551 *985 *364 49 827 *384 85 265 587 421]
% 
% ⭐Bin Cell -2-
% [810 *1008 386 614 850 *849 126 1022 367 677 347 284 27 551 991 990 49 827 822 85 265 587 421]

toc;
fclose all;
