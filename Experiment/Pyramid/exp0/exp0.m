% 画出所选文件的STFT，并给出bin值
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

[loraSet] = readLoraSet('sf8_BW125.json');
% 读取原始信号文件
readDir = strcat('E:\Pyramid_samples\Samples_base\BW', string(loraSet.bw), '\sf', string(loraSet.sf), '\');
% 合成冲突文件位置
writeDir = 'E:\Pyramid_samples\conflict\';
% bin参考文件夹位置
refDir = 'E:\Pyramid_samples\';
% 设置保存结果的文件名
saveName = strcat('exp0_sf', string(loraSet.sf), '_bw', string(loraSet.bw), '.mat');

fileIn = dir(fullfile(readDir, '*.sigmf-data')); 
% 节点数量
nodeNumArray = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000];

flagMatlabName = strcat(writeDir, 'flagMatlab.txt');
flagGrloraName = strcat(writeDir, 'flagGrlora.txt');
resultArray = zeros(1, length(nodeNumArray));

% 获取节点信息
load(strcat(refDir, 'node_info.mat'));

for nodeIndex = 1:length(nodeNumArray)
    resultAll = 0;
    numAll = 0;
    srrBefore = 0;
    nodeNum = nodeNumArray(nodeIndex);
    nodeInfoTmp = node_info(1:nodeNum, :);
    nodeSelect = find(nodeInfoTmp(:,1) == 1);   % 选择在1信道的节点
    nodeInfoTmp = nodeInfoTmp(nodeSelect, :);
    while true
        % 根据参数生成发包时间矩阵
        pkgRate = 4;
        pktNum = 5;
        pktLength = loraSet.dine * (12.25+40);
        fs = loraSet.sample_rate;
        [T] = gen_pkt_time_pyramid(nodeNum, pkgRate, pktNum, pktLength, fs);
        % 根据节点信息挑选节点
        T = T(nodeSelect, :);
        T = T(:, 2:end);
        lengthAll = max(T,[],"all");
        samples = zeros(1, ceil(lengthAll) + (40+12.25)*loraSet.dine);
        
        % 合成冲突信号
        for row = 1:size(T,1)
            fail = nodeInfoTmp(row, end);
            for col = 1:size(T,2)
                % 随机取信号
                id = ceil(length(fileIn)*rand());
                signal = readSignalFile(readDir, fileIn(id));
                % 衰落信号
                signal = sqrt(10^(fail/10)).*signal;
                % 获得时间戳
                timeStamp = ceil(T(row, col));
                % 叠加信号
                samplesTmp = samples(timeStamp+1 : timeStamp+(40+12.25)*loraSet.dine);
                samplesTmp = samplesTmp + signal(1:(40+12.25)*loraSet.dine);
                samples(timeStamp+1 : timeStamp+(40+12.25)*loraSet.dine) = samplesTmp;
            end
        end
        samples = awgn(samples, 20, 'measured');
        
        % 写入文件中
        signalLength = length(samples)*2;
        signalProcessed = zeros(1, signalLength);
        signalReal = real(samples);   signalImag = imag(samples);
        signalProcessed(1:2:signalLength-1) = signalReal;
        signalProcessed(2:2:signalLength) = signalImag;
        
        writePath = strcat(writeDir, 'conflict.sigmf-data');
        fid=fopen(writePath, 'wb');
        fwrite(fid, signalProcessed, 'float32');
        fclose all;
        % 写入标志位文件
        writematrix([],flagMatlabName);

        fprintf("wait gr-lora...\n");
        % 等待gr-lora处理完成
        while ~exist(flagGrloraName,'file')
        end
        % 清除gr-lora标志文件
        delete(flagGrloraName);

        % gr-lora处理完成，统计结果
        fprintf("gr-lora Done!\n");
        bin_ref = load(strcat(refDir, 'bin_ref', string(loraSet.sf), '.txt'))';
        bin_file = strcat(refDir, 'bin.txt');
        if exist(bin_file, 'file')
            bin = textread(bin_file,'','delimiter',',','emptyvalue',NaN);
        else
            bin = 10000*ones(1, length(bin_ref));
        end
        fclose all;
        
        % 检查bin的长度是否满足要求
        if size(bin,2) < length(bin_ref)
            bin = [bin(), 10000*ones(size(bin), length(bin_ref) - size(bin,2))];
        end
        pkg_num = length(T(:));
        % 统计结果
        for i = 1:size(bin, 1)
            bin_tmp = bin(i, 1:length(bin_ref));
            tmp = (bin_tmp - (bin_ref-1));
            resultAll = resultAll + sum((bin_tmp - (bin_ref-1)) == 0);
        end
        numAll = numAll + pkg_num * length(bin_ref);
        % 计算当前实验总SRR
        srrNow = resultAll/numAll;
        fprintf("node num is %d srrNow is %f\n", nodeNumArray(nodeIndex), srrNow);
        % 如果收敛，退出循环
        if abs((srrNow - srrBefore)/srrBefore) <= 0.01
            break;
        else
            srrBefore = srrNow;
        end
        
    end
    resultArray(nodeIndex) = srrNow;
    save(saveName, "resultArray");
end

toc;
fclose all;