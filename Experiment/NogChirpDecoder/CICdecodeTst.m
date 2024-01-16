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
fileDir = '\\192.168.3.102\e\data\ChNum_2_m2h3\';
fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
% 从文件中读取信号流
[signal] = readSignalFile(fileDir, fileIn(1));

CICDecoder = CICDecoder.decode(signal);
for i = 1 : 1: length(CICDecoder.binRecord)
    disp(CICDecoder.binRecord(i));
end


% [810 722 386 614 850 593 126 1022 842 677 347 862 581 551 574 113 806 827 822 85 421 587 421]
% [810 1010 386 614 850 406 126 1022 367 677 347 284 27 551 1004 1006 49 836 822 85 265 587 421]
% 准确率 14/23 20/23  0/23

toc;
fclose all;
