% 测试CHchirpMultiDecoder解单包的性能，保存实验结果
clear;
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 所有采样信号的基本参数信息枚举
bwArr = [125e3, 250e3];
sfArr = [9, 10];
channelNumArr = [2, 4, 6];
subchirpNumArr = [2, 4];
samplesRate = 2e6;
DEBUG = false;

% 遍历所有配置，执行实验
for bw = bwArr
    for sf = sfArr
        for channelNum = channelNumArr
            for subchirpNum = subchirpNumArr
                % 读取配置和验证文件
                [loraSet] = readLoraSet('GeneralConfig.json', sf, bw, samplesRate);
                loraSet.channelNum = channelNum; % 信道数目
                loraSet.subchirpNum = subchirpNum; % subchirp数目
                loraSet.payloadNum = 20; % payload数目
                clear CHchirpMultiDecoder;
                CHchirpMultiDecoder = CHchirpMultiDecoder(loraSet, DEBUG);
                % 从文件中加载bin groundtrurh
                true_bin = importdata(strcat('.\Code\Config\bin\SF', string(sf), '.txt'))';
                % 读取文件夹下所有采样值文件
                fileDir = 'D:\CHchirp23_11\SF' + string(loraSet.sf) ...
                    + '_BW' + string(loraSet.bw/1000) ...
                    + 'K_sub' + string(loraSet.subchirpNum) ...
                    + '_channels' + string(loraSet.channelNum);
                if exist(fileDir, "dir")
                    record = cell(3, 1);  % 设置元胞数组记录结果
                    argsName = 'SF' + string(loraSet.sf) ...
                        + '_BW' + string(loraSet.bw/1000) ...
                        + '_subchirp' + string(loraSet.subchirpNum) ...
                        + '_channel' + string(loraSet.channelNum);
                    record{1} = argsName;
                    fileDir = fileDir + "\";
                    fileIn = dir(fullfile(fileDir, '*.sigmf-data'));
                    fileLength = length(fileIn);
                    recordBin = cell(1, fileLength);
                    recordRate = zeros(1, fileLength);
                    % 读取该文件夹下所有的文件
                    for fileCount = 1:fileLength
                        % 从文件中读取信号流
                        [signal] = readSignalFile(fileDir, fileIn(fileCount));
                        try
                            CHchirpMultiDecoder = CHchirpMultiDecoder.decode(signal);
                        catch
                            disp(argsName + "，文件编号fileCount：" + fileCount + "报错");
                        end
                        % 验证解调出来的所有bin是否正确
                        trueRateMax = 0;
                        for binResultIndex = 1:length(CHchirpMultiDecoder.binArray)
                            [true_chirp, true_rate] = vertify_bin(CHchirpMultiDecoder.binArray{binResultIndex}, true_bin);
                            trueRateMax = max(trueRateMax, true_rate);
                            if trueRateMax == 1
                                recordBin{fileCount} = CHchirpMultiDecoder.binArray{binResultIndex};
                                break;
                            end
                            recordBin{fileCount} = CHchirpMultiDecoder.binArray{binResultIndex};
                        end
                        disp(argsName + " 文件序号:" + fileCount + ", 正确率为：" + trueRateMax*100 + "%");
                        recordRate(fileCount) = trueRateMax;
                    end
                    % 记录结果
                    record{2} = recordBin;
                    record{3} = recordRate;
                    % 画出直方图，并保存直方图
                    binEdges = 0:0.01:1;  % 这里设置了20个柱子，可以根据需要调整
                    histogram(recordRate, binEdges, 'Normalization', 'probability');
                    title(argsName+'成功率分布');
                    xlabel('成功率');
                    ylabel('概率');
                    saveas(gcf, '.\Result\singlePackageExp\' + argsName + '.fig');
                    saveas(gcf, '.\Result\singlePackageExp\' + argsName + '.png');
                    save('.\Result\singlePackageExp\' + argsName + '.mat', 'record');
                end
            end
        end
    end
end


toc;
fclose all;