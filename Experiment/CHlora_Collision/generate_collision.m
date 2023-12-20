% 画出所选文件的STFT，并给出bin值
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 采样值文件读取路径和保存路径
inDir = 'E:\share\aligned\';
% 读取配置和验证文件
[loraSet] = readLoraSet('sf10_BW125.json', 2e6);
% 读取文件夹下所有采样值文件
fileIn = dir(fullfile(inDir, '*.sigmf-data'));
signal_all = zeros(1, 120*loraSet.dine);
len = length(signal_all);
% for fileCount = 1:length(fileIn)
for fileCount = 1:3
    % 从文件中读取信号流
    [signal] = readSignalFile(inDir, fileIn(fileCount));
    offset = round(rand(1)*10*loraSet.dine);
    signal = [zeros(1, offset), signal, zeros(1, 120*loraSet.dine-length(signal)-offset)];
    signal_all = signal_all + signal;
end
write_signal_to_file(signal_all, strcat('E:\share\collision\', 'collision', '.sigmf-data'));
toc;
fclose all;     %关闭所有matlab打开的文件