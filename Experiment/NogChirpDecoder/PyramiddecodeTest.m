% 
clear;
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 所有采样信号的基本参数信息枚举
bw = 125e3;
sf = 10;
samplesRate = 2e6;

% 读取配置和验证文件
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
loraSet.payloadNum = 23; % payload数目
SignalLength = loraSet.dine*80;  % 整个信号的最大长度

% 初始化decoder
PyramidDecoder = PyramidDecoder(loraSet);

% fileDir = 'd:\data\ChNum_2_m2h3\';
% fileDir = 'd:\data\Collision-2_CH-2\';

fileDir = 'D:\data\SameBinInterfer\';
fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
[signal1] = readSignalFile(fileDir, fileIn(2));
fileDir = 'D:\data\SameBinInterfer\';
fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
[signal2] = readSignalFile(fileDir, fileIn(1));

signal1 = [signal1 zeros(1, 50000)];
signal2 = [zeros(1, 50000) signal2];

signal = signal1 + signal2;


emptySignal = zeros(1, 10000); % create an array of zeros with the specified length
paddedSignal = [emptySignal signal]; % concatenate the empty signal with the original signal

PyramidDecoder = PyramidDecoder.decode(paddedSignal);

dimensions = size(PyramidDecoder.binRecord);
lengths = cellfun(@length, PyramidDecoder.binRecord);
disp("⭐payloadBin dimensions: " + num2str(dimensions(1))+ "x" + num2str(dimensions(2)));
disp("⭐payloadBin lengths: " + num2str(lengths));

for i = 1:numel(PyramidDecoder.binRecord)
    fprintf('\n⭐Bin Cell -%d-\n', i);
    % 使用 cellfun 将当前行的每个元素格式化为字符串，并连接起来
    row_str = cellfun(@(x) sprintf('%s', mat2str(x)), PyramidDecoder.binRecord(i), 'UniformOutput', false);
    % 使用 strjoin 将格式化后的字符串连接起来，并输出
    disp(strjoin(row_str, ' '));

end

toc;
fclose all;
