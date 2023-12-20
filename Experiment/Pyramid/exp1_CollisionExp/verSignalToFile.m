% 测试CHchirpDecoder
clear;
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 读取配置和验证文件
sf = 9;
bw = 125e3;
samplesRate = 2e6;
% 设置计数器
SuccessNum = 1;
[loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
loraSet.channelNum = 2; % 信道数目
loraSet.subchirpNum = 2; % subchirp数目
loraSet.payloadNum = 20; % payload数目
PyramidDecoder = PyramidDecoder(loraSet);
% 读取文件夹下所有采样值文件
fileDir = "D:\CHchirp_20Nodes\tmp\";
writeDir = "D:\CHchirp_20Nodes\tmp\write\";
fileDir = fileDir + 'SF' + string(sf) + 'BW' + (bw/1000) + '\';
writeDir = writeDir + 'SF' + string(sf) + 'BW' + (bw/1000) + '\';
% 读取正确的bin值
true_bin = importdata(strcat('.\Code\Config\bin\SF', string(sf), '.txt'))';


for RLTIndex = 1:5
    filePath = fileDir + "RTL" + string(RLTIndex) + "\";
    fileIn = dir(fullfile(filePath, '*.sigmf-data'));
    count = 0;
    for fileCount = 1:length(fileIn)
        % 从文件中读取信号流
        trueRate = 0;
        [signal] = readSignalFile(filePath, fileIn(fileCount));
        signalTmp = [zeros(1, 2*loraSet.dine), signal, zeros(1, 10*loraSet.dine)];
    
        % 写入读取区
        try
            PyramidDecoder = PyramidDecoder.decode(signalTmp);
        catch
            continue;
        end
        % 计算正确率
        if ~isempty(PyramidDecoder.binRecord)
            trueRate = calculateAccuracy(true_bin, PyramidDecoder.binRecord{1});
        end
        
        disp("当前文件序号：" + fileCount);
        disp(PyramidDecoder.binRecord);
        disp("正确率: " + trueRate*100 + "%");
        if trueRate == 1
            count = count + 1;
            fileName = writeDir + "RTL" + string(RLTIndex) + "\" + string(count) + ".sigmf-data";
            write_signal_to_file(signal, fileName);
            if count == SuccessNum
                break;
            end
        end
    end
end

toc;
fclose all;