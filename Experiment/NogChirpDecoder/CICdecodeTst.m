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
% fileDir = 'd:\data\ChNum_2_m2h3\';
% fileDir = 'd:\data\Collision-2_CH-2\';
fileDir = '\\192.168.3.102\e\data\ChNum_2_m2h3_22\';
% fileDir = 'd:\data\ChNum_2_m2h3\';
fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
true_bin = importdata(strcat('.\Config\bin\NogSF', string(sf), '.txt'))';

% 从文件中读取信号流
bar = waitbar(0, 'Loading...');    % waitbar 显示进度条

fileNumber = numel(fileIn);
fileBinTrueRate = zeros(1, fileNumber);
for file_i = 1 : fileNumber
    [signal] = readSignalFile(fileDir, fileIn(file_i));
    emptySignal = zeros(1, 10000); % create an array of zeros with the specified length
    paddedSignal = [emptySignal signal]; % concatenate the empty signal with the original signal

    try
        CICDecoder = CICDecoder.decode(paddedSignal);
    catch
        disp("文件: " + num2str(file_i) + " 出现错误");
    end
    % 计算该文件的正确率
    pktNum = length(CICDecoder.binRecord);
    if pktNum
        % 计算该文件的正确率
        recordRateTmp = zeros(1, pktNum);
        for binResultIndex = 1 : pktNum
            recordRateTmp(binResultIndex) = calculateAccuracy(true_bin, CICDecoder.binRecord{binResultIndex});
        end
        recordRateTmp = sort(recordRateTmp, 'descend');
        fileBinTrueRate(file_i) = mean(recordRateTmp);
        disp(['文件 ', num2str(file_i) , ' 检测到 ', num2str(pktNum), ' 个信号', ' •准确率: ', num2str(fileBinTrueRate(file_i) * 100), '%']);
    else
        disp(['文件 ', num2str(file_i) , ' 未检测到信号']);
    end
    waitbarStr = ['目前进度 ', num2str(100 * file_i / fileNumber), '%, 已完成 ', num2str(file_i), '/', num2str(fileNumber)];   % 显示的文本
    waitbar(file_i / fileNumber, bar, waitbarStr);
end
close(bar);

fileBinTrueRate = fileBinTrueRate(fileBinTrueRate ~= 0);
disp(['综合准确率: ', num2str(mean(fileBinTrueRate)*100), '%']);
% 综合准确率: 90.7971%
% 历时 260.721951 秒。


% emptySignal = zeros(1, 10000); % create an array of zeros with the specified length
% paddedSignal = [emptySignal signal]; % concatenate the empty signal with the original signal
% 
% CICDecoder = CICDecoder.decode(paddedSignal);
% 
% dimensions = size(CICDecoder.binRecord);
% lengths = cellfun(@length, CICDecoder.binRecord);
% disp("⭐payloadBin dimensions: " + num2str(dimensions(1))+ "x" + num2str(dimensions(2)));
% disp("⭐payloadBin lengths: " + num2str(lengths));
% 
% for i = 1:numel(CICDecoder.binRecord)
%     fprintf('\n⭐Bin Cell -%d-\n', i);
%     % 使用 cellfun 将当前行的每个元素格式化为字符串，并连接起来
%     row_str = cellfun(@(x) sprintf('%s', mat2str(x)), CICDecoder.binRecord(i), 'UniformOutput', false);
%     % 使用 strjoin 将格式化后的字符串连接起来，并输出
%     disp(strjoin(row_str, ' '));
%     row_str = replace(row_str, '[', '');
%     row_str = replace(row_str, ']', '');
%     % 计算准确率
%     accuracy = sum(str2num(row_str{1}) == true_bin) / numel(true_bin) * 100; % 将匹配的数量除以总数，并乘以100以获得百分比
%     disp(['准确率：', num2str(accuracy), '%']);
% end

% accuracy = true_bin_Num / binSum * 100;
% disp("正确解出的 Bin 值: " + num2str(true_bin_Num));
% disp("Bin 值总数: " + num2str(binSum));
% disp(['准确率: ', num2str(accuracy), '%']);

% 正确解出的 Bin 值: 735
% Bin 值总数: 805
% 准确率: 91.3043%
% 历时 196.450467 秒。

% ⭐Bin Cell -1-
% [810 1010 386 614 850 406 126 1022 367 677 347 284 27 551 *985 *364 49 827 *384 85 265 587 421]
% 
% ⭐Bin Cell -2-
% [810 *1008 386 614 850 *849 126 1022 367 677 347 284 27 551 991 990 49 827 822 85 265 587 421]

toc;
fclose all;
