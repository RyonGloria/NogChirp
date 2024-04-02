% 测试重叠信号解码
clear;
fclose all;     % 关闭所有 matlab 打开的文件
tic;            % 打开计时器

% 基本参数设置
sf = 10;
bw = 250e3;
samplesRate = 2e6;
debugPath = "terminal";
DebugLevel = 4;
DebugUtil = DebugUtil(DebugLevel, debugPath);

% 读取配置和验证文件`
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
loraSet.payloadNum = 23;  % payload数目

%% 读取文件夹下所有采样值文件
% fileDir = '\\192.168.3.102\e\data\ChNum_2_m2h3\';
% fileDir = 'd:\data\ChNum_2_m2h3\';
% fileDir = '\\192.168.3.102\e\share\samples\';
% fileDir = 'd:\data\Collision-2_CH-2\';
fileDir = 'd:\data\1_17indoor\FFT_jun\';
% fileDir = 'd:\data\ChNum_2_m2h3_22\';
fileIn = dir(fullfile(fileDir, '*.sigmf-data'));

% true_bin = importdata(strcat('.\Config\bin\NogSF', string(sf), '.txt'))';


% 从文件中读取信号流
fileNum = numel(fileIn);
fileBinTrueRate = zeros(1, fileNum);
for file_i = 1 : fileNum
    [signal] = readSignalFile(fileDir, fileIn(file_i));
    obj = VarCutChirpDecoder(loraSet, DebugUtil);
    try
        obj = obj.decode(signal);
    catch
        disp(['文件 ' , num2str(file_i), '/', num2str(fileNum), ' 出现错误']);
        continue;
    end

    pktNum = length(obj.payloadBin);
    if pktNum
        for i = 1:numel(obj.payloadBin)
        row_str = cellfun(@(x) sprintf('%s', mat2str(x)), obj.payloadBin{i}, 'UniformOutput', false);
        disp(strjoin(row_str, ' '));
        end
    else
        disp(['文件 ', num2str(file_i) , '/', num2str(fileNum), ' 未检测到信号']);
    end
    % if pktNum
    %     % 计算该文件的正确率
    %     recordRateTmp = zeros(1, pktNum);
    %     for binResultIndex = 1 : pktNum
    %         recordRateTmp(binResultIndex) = calculateAccuracy(true_bin, cell2mat(obj.payloadBin{binResultIndex}));
    %     end
    %     recordRateTmp = sort(recordRateTmp, 'descend');
    %     fileBinTrueRate(file_i) = mean(recordRateTmp);
    %     disp(['文件 ', num2str(file_i) , '/', num2str(fileNum), ' 检测到 ', num2str(pktNum), ' 个信号', ' (准确率: ', num2str(fileBinTrueRate(file_i) * 100), '%)']);
    % else
    %     disp(['文件 ', num2str(file_i) , '/', num2str(fileNum), ' 未检测到信号']);
    % end
end

% fileBinTrueRate = fileBinTrueRate(fileBinTrueRate ~= 0);
% disp(['综合准确率: ', num2str(mean(fileBinTrueRate)*100), '%']);

% 综合准确率: 97.0109%
% 历时 19.772260 秒。

% for file_i = 1 : fileNum
    % [signal] = readSignalFile(fileDir, fileIn(file_i));
    % obj = VarCutChirpDecoder(loraSet, DebugUtil);
    % obj = obj.decode(signal);
    % dimensions = size(obj.payloadBin);
    % Get the length of each dimension
    % lengths = cellfun(@length, obj.payloadBin);
    % disp("⭐payloadBin dimensions: " + num2str(dimensions(1))+ "x" + num2str(dimensions(2)));
    % disp("⭐payloadBin lengths: " + num2str(lengths));
    % fprintf('\n');
    % 循环遍历每个单元格
    % for i = 1:numel(obj.payloadBin)
        % disp("⭐Bin Cell: Channel -" + obj.detectedPktAll{i}(1) + "-, PreambleEndPos -" + obj.detectedPktAll{i}(2) + "-")
        % 使用 cellfun 将当前行的每个元素格式化为字符串，并连接起来
        % row_str = cellfun(@(x) sprintf('%s', mat2str(x)), obj.payloadBin{i}, 'UniformOutput', false);
        % 使用 strjoin 将格式化后的字符串连接起来，并输出
        % disp(strjoin(row_str, ' '));
        % accuracy = sum(str2double(row_str) == true_bin) / numel(true_bin) * 100; % 将匹配的数量除以总数，并乘以100以获得百分比
        % disp(['准确率：', num2str(accuracy), '%']);
        % fprintf('\n');
    % end
    % waitbarStr = ['目前进度 ', num2str(100 * file_i / fileNum), '%, 已完成 ', num2str(file_i), '/', num2str(fileNum)];   % 显示的文本
    % waitbar(file_i/fileNum, bar, waitbarStr);
    % fprintf('\n');
% end

% 正确解出的 Bin 值: 1051
% Bin 值总数: 1081
% 准确率: 97.2248%
% 历时 20.209275 秒。
toc;
fclose all;
