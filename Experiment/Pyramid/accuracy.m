% 测试CHchirpDecoder
clear;
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 读取配置和验证文件
sf = 10;
bw = 125e3;
samplesRate = 2e6;
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
loraSet.channelNum = 2; % 信道数目
loraSet.subchirpNum = 2; % subchirp数目
loraSet.payloadNum = 20; % payload数目
PyramidDecoder = PyramidDecoder(loraSet);
% 读取文件夹下所有采样值文件
fileDir = 'E:\share\samples\';
% fileDir = 'D:\CHchirp_20Nodes\Pyramid_SF10_BW125K\Node1\tmp\';
writeDir = 'D:\CHchirp_20Nodes\USRP\Pyramid_SF10_BW125K\RTL1\';
true_bin = importdata(strcat('.\Code\Config\bin\SF', string(sf), '.txt'))';

fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
for fileCount = 1:length(fileIn)
% for fileCount = 1:1
    % 从文件中读取信号流
    [signal] = readSignalFile(fileDir, fileIn(fileCount));
    signalTmp = [zeros(1, loraSet.dine*10), signal, zeros(1, loraSet.dine*40)];
    % 写入读取区
    PyramidDecoder = PyramidDecoder.decode(signalTmp);
    
    disp("当前文件序号：" + fileCount);
    disp(PyramidDecoder.binRecord);
    if ~isempty(PyramidDecoder.binRecord)
        trueRate = calculateAccuracy(true_bin, PyramidDecoder.binRecord{1});
        disp("正确率: " + trueRate*100 + "%");
        fileName = writeDir + string(count) + ".sigmf-data";
        write_signal_to_file(signal, fileName);
    end
end

toc;
fclose all;