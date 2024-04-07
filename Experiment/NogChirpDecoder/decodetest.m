% 测试重叠信号解码
clear;
fclose all;     % 关闭所有 matlab 打开的文件
tic;            % 打开计时器

% 基本参数设置
sf = 10;
bw = 125e3;
samplesRate = 2e6;

%% 读取配置和验证文件
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
loraSet.payloadNum = 23;  % payload数目
%% Decoder
obj = NogChirpDecoder(loraSet);

%% Symbols convert to Bytes
phy = LoRaPhyLay(sf, bw, samplesRate);
phy.has_header = 0;         % explicit header mode
phy.cr = 1;                 % code rate = 4/8 (1:4/5 2:4/6 3:4/7 4:4/8)
phy.crc = 1;                % enable payload CRC checksum
phy.payload_len = 17;       % payload length 17 bytes
phy.is_debug = 0;           % enable debug mode

%% 读取文件夹下所有采样值文件
% fileDir = '\\192.168.3.102\e\data\ChNum_2_m2h3\';
% fileDir = 'd:\data\ChNum_2_m2h3\';
% fileDir = '\\192.168.3.102\e\share\samples\';
fileDir = 'd:\data\ChNum_1_m2\';

fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
true_bin = importdata(strcat('.\Config\bin\NogSF', string(sf), '.txt'))';
% 从文件中读取信号流
[signal] = readSignalFile(fileDir, fileIn(2));

%% Decode Two Channel
obj = obj.decodeTwoCH(signal);
dimensions = size(obj.payloadBin);
% Get the length of each dimension
lengths = cellfun(@length, obj.payloadBin);
disp("⭐payloadBin dimensions: " + num2str(dimensions(1))+ "x" + num2str(dimensions(2)));
disp("⭐payloadBin lengths: " + num2str(lengths));
fprintf('\n');

BinRes = zeros(1, obj.loraSet.payloadNum);
% 循环遍历每个单元格
for i = 1:numel(obj.payloadBin)
    disp("⭐Bin Cell: Channel -" + obj.ResultInfo{i}(1) + "-, PreambleStartPos -" + obj.ResultInfo{i}(2) + "-")
    % 使用 cellfun 将当前行的每个元素格式化为字符串，并连接起来
    row_str = cellfun(@(x) sprintf('%s', mat2str(x)), obj.payloadBin{i}, 'UniformOutput', false);
    % 使用 strjoin 将格式化后的字符串连接起来，并输出
    disp(strjoin(row_str, ' '));
    fprintf("\n");
    [data, checksum] = phy.decode(cell2mat(obj.payloadBin{i})');
    fprintf("[decode] bytes:\n");
    disp(data(1:17)');
    fprintf("[decode] checksum:\n");
    disp(checksum');
    data_str = char(data(1:17)');
    fprintf("[decode] string:\n");
    disp(data_str);
    % 计算准确率
    accuracy = sum(str2double(row_str) == true_bin) / numel(true_bin) * 100; % 将匹配的数量除以总数，并乘以100以获得百分比
    % disp(['准确率：', num2str(accuracy), '%']);
    fprintf('\n');
end
fprintf('\n');

toc;
fclose all;

% ======================================preamble bin: [980]======================================
% preambleEndPos: 9
% ======================================preamble bin: [958]======================================
% preambleEndPos: 21