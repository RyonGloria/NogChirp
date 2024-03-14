% 测试重叠信号解码
clear;
fclose all;     % 关闭所有 matlab 打开的文件
tic;            % 打开计时器

% 基本参数设置
sf = 10;
bw = 125e3;
samplesRate = 2e6;

% 读取配置和验证文件
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
loraSet.payloadNum = 23;  % payload数目
obj = NogChirpDecoder(loraSet);
%% 读取文件夹下所有采样值文件
% fileDir = '\\192.168.3.102\e\data\ChNum_2_m2h3\';
fileDir = 'd:\data\ChNum_2_m2h3\';
fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
% 从文件中读取信号流
[signal] = readSignalFile(fileDir, fileIn(1));

%% Decode Two Channel
obj = obj.decodeTwoCH(signal);
dimensions = size(obj.payloadBin);
% Get the length of each dimension
lengths = cellfun(@length, obj.payloadBin);
disp("⭐payloadBin dimensions: " + num2str(dimensions(1))+ "x" + num2str(dimensions(2)));
disp("⭐payloadBin lengths: " + num2str(lengths));
% 循环遍历每个单元格
for i = 1:numel(obj.payloadBin)
    fprintf('\n');
    disp(['⭐Bin Cell -', num2str(i), '-']);
    fprintf('%5d', obj.payloadBin{i}); % 控制每个数字的宽度为 10
end
fprintf('\n');


toc;
fclose all;

% ======================================preamble bin: [980]======================================
% preambleEndPos: 9
% ======================================preamble bin: [958]======================================
% preambleEndPos: 21