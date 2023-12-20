classdef CHchirpMultiDecoder < CHchirpDecoder
    properties
        % 子解调器
        CHchirpCollisionDecoder;
        % 检测队列
        detect_array;           % 存放4个信道检测到的bin
        detect_array_count;     % 计数器
        detect_array_number;    % 记录四个信道队列中存放有效数据的数目
        % 结果集
        binArray;
        detectCount;
        detectSFDCount;
        SFDPeakAmp;
        % Debug
        DebugUtil;
    end
    
    methods
        function obj = CHchirpMultiDecoder(loraSet, DebugUtil)
            obj@CHchirpDecoder(loraSet);
            obj.DebugUtil = DebugUtil;
            obj.CHchirpCollisionDecoder = CHchirpCollisionDecoder(loraSet, obj.channelMatrix, DebugUtil);
            % 检测队列
            obj.detect_array = zeros(obj.channelNum, loraSet.fft_x);   % 存放每个信道检测到的峰值bin（bin最多为fft_x种）
            obj.detect_array_count = zeros(obj.channelNum, loraSet.fft_x);  % 计数器
            obj.detect_array_number = zeros(1, obj.channelNum);  % 记录四个信道队列中存放有效数据的数目
            % 结果集
            obj.binArray = cell(0, 1);
        end

        function obj = clearQueue(obj)
            % 检测队列
            obj.detect_array = zeros(obj.channelNum, obj.loraSet.fft_x);   % 存放每个信道检测到的峰值bin（bin最多为fft_x种）
            obj.detect_array_count = zeros(obj.channelNum, obj.loraSet.fft_x);  % 计数器
            obj.detect_array_number = zeros(1, obj.channelNum);  % 记录四个信道队列中存放有效数据的数目
            % 结果集
            obj.binArray = cell(0, 1);
            obj.SFDPeakAmp = cell(0,1);
            obj.detectCount = 0;
            obj.detectSFDCount = 0; 
        end

        function obj = decode(obj, signals)
            obj = obj.clearQueue();
            % 每次调用decode解码信号前，需要将检测队列清空
            dine = obj.loraSet.dine; % 每个chirp的长度
            windowsNum = fix(length(signals) / dine); % 计算需要移动多少个窗口进行信号检测和解调

            for window = 1:windowsNum % 开始扫描每一个窗口，发现其中的premable
                %{
                    划分出每个信道的信号（这里重写了CHchirpDecoder的divideChannel函数，因为要实现多解码无法将划分信道的信号记录在类变量中
                        所以以输入输出信号的方式设计）
                %}
                windowSignal = signals((window-1)*dine+1 : window*dine);
                splitSignal = obj.divideChannel(windowSignal);
                % 对每个信道进行信号检测
                for channel = 1:obj.channelNum
                    % 对每个信道进行检测，是否存在lora数据包，并标记preamble对应的bin
                    [obj, preambleBin] = obj.detect(channel, splitSignal);
                    % 存在满足7个窗口峰值bin相同的preambleBin
                     if ~isempty(preambleBin)
                        % 根据detect的结果进行解调
                        for num = 1:length(preambleBin)
                            % disp(preambleBin(num));
                            % 重新初始化CHchirpCollisionDecoder
                            obj.CHchirpCollisionDecoder = obj.CHchirpCollisionDecoder.setArgs(channel, window, preambleBin(num));
                            obj.CHchirpCollisionDecoder = obj.CHchirpCollisionDecoder.decode(signals);
                            if obj.CHchirpCollisionDecoder.errorFlag == false
                                obj.detectSFDCount = obj.detectSFDCount + 1;
                                obj.SFDPeakAmp{obj.detectSFDCount} = obj.CHchirpCollisionDecoder.SFDPeakAmp(1) / obj.CHchirpCollisionDecoder.SFDPeakAmp(2);
                                obj.binArray{end+1} = obj.CHchirpCollisionDecoder.payloadBin;
                            end
                            obj.detectCount = obj.detectCount + 1;
                        end
                    end
                end
            end
        end

        function obj = decodeCalculateDetect(obj, signals)
            obj = obj.clearQueue();
            % 每次调用decode解码信号前，需要将检测队列清空
            dine = obj.loraSet.dine; % 每个chirp的长度
            windowsNum = fix(length(signals) / dine); % 计算需要移动多少个窗口进行信号检测和解调

            % 临时变量，记录数据
            for window = 1:windowsNum % 开始扫描每一个窗口，发现其中的premable
                %{
                    划分出每个信道的信号（这里重写了CHchirpDecoder的divideChannel函数，因为要实现多解码无法将划分信道的信号记录在类变量中
                        所以以输入输出信号的方式设计）
                %}
                windowSignal = signals((window-1)*dine+1 : window*dine);
                splitSignal = obj.divideChannel(windowSignal);
                % 对每个信道进行信号检测
                for channel = 1:obj.channelNum
                    % 对每个信道进行检测，是否存在lora数据包，并标记preamble对应的bin
                    [obj, preambleBin] = obj.detect(channel, splitSignal);
                    % 存在满足7个窗口峰值bin相同的preambleBin
                    if ~isempty(preambleBin)
                        % 根据detect的结果进行解调
                        for num = 1:length(preambleBin)
                            % disp(preambleBin(num));
                            % 重新初始化CHchirpCollisionDecoder
                            obj.CHchirpCollisionDecoder = obj.CHchirpCollisionDecoder.setArgs(channel, window, preambleBin(num));

                            % obj.CHchirpCollisionDecoder = obj.CHchirpCollisionDecoder.decode(signals);
                            obj.CHchirpCollisionDecoder = obj.CHchirpCollisionDecoder.detectPass(signals);
                            if obj.CHchirpCollisionDecoder.errorFlag == false
                                obj.detectSFDCount = obj.detectSFDCount + 1;
                                obj.SFDPeakAmp{obj.detectSFDCount} = obj.CHchirpCollisionDecoder.SFDPeakAmp(1) / obj.CHchirpCollisionDecoder.SFDPeakAmp(2);
                            end
                            obj.detectCount = obj.detectCount + 1;
                        end
                        
                    end
                end
            end
        end

        function obj = decodeDebug(obj, signals)
            obj = obj.clearQueue();
            % 每次调用decode解码信号前，需要将检测队列清空
            dine = obj.loraSet.dine; % 每个chirp的长度
            windowsNum = fix(length(signals) / dine); % 计算需要移动多少个窗口进行信号检测和解调

            % 临时变量，记录数据
            for window = 1:windowsNum % 开始扫描每一个窗口，发现其中的premable
                %{
                    划分出每个信道的信号（这里重写了CHchirpDecoder的divideChannel函数，因为要实现多解码无法将划分信道的信号记录在类变量中
                        所以以输入输出信号的方式设计）
                %}
                windowSignal = signals((window-1)*dine+1 : window*dine);
                splitSignal = obj.divideChannel(windowSignal);
                % 对每个信道进行信号检测
                for channel = 1:obj.channelNum
                    % 对每个信道进行检测，是否存在lora数据包，并标记preamble对应的bin
                    [obj, preambleBin] = obj.detect(channel, splitSignal);
                    % 存在满足7个窗口峰值bin相同的preambleBin
                    if ~isempty(preambleBin)
                        % 根据detect的结果进行解调
                        for num = 1:length(preambleBin)
                            % disp(preambleBin(num));
                            % 重新初始化CHchirpCollisionDecoder
                            obj.CHchirpCollisionDecoder = obj.CHchirpCollisionDecoder.setArgs(channel, window, preambleBin(num));

                            % obj.CHchirpCollisionDecoder = obj.CHchirpCollisionDecoder.decode(signals);
                            obj.CHchirpCollisionDecoder = obj.CHchirpCollisionDecoder.decodeDebug(signals);
                            if obj.CHchirpCollisionDecoder.errorFlag == false
                                obj.detectSFDCount = obj.detectSFDCount + 1;
                                obj.SFDPeakAmp{obj.detectSFDCount} = obj.CHchirpCollisionDecoder.SFDPeakAmp(1) / obj.CHchirpCollisionDecoder.SFDPeakAmp(2);
                                obj.binArray{end+1} = obj.CHchirpCollisionDecoder.payloadBin;
                            end
                            obj.detectCount = obj.detectCount + 1;
                        end
                        
                    end
                end
            end
        end

        %{
            检测信道中是否存在满足preamble特征的信号
            channel：信道的标号
            splitSignal：划分信道后的信号
            preambleBin：满足特征的信号preambleBin集合
        %}
        function [obj, preambleBin] = detect(obj, channel, splitSignal)
            preambleBin = [];  % 用于存放检测出来的lora包的preamble bin
            dine = obj.loraSet.dine;
            fft_x = obj.loraSet.fft_x;
            leakWidth = obj.loraSet.leakage_width1;
            if isActive(obj, splitSignal(channel, :))  % 如果这个信道存在信号
                % 对这个信道的信号乘以downchirp，做FFT，通过findpeaks函数找到突出峰值对应的bin，来判断preamble特性
                signal_tmp = splitSignal(channel, :);
                dechirp = signal_tmp .* obj.idealDownchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine); 
                % 找出若干个峰值，每个峰值间隔需超过leakWidth
                result = obj.findpeaksWithShift(dechirp_fft, fft_x);
                binPos = result(2, :);
                % [~, binPos] = findpeaks(dechirp_fft, "MinPeakDistance", fft_x*leakWidth, "MinPeakHeight", mean(dechirp_fft) * 10);
                % 对每一个筛选过的峰值进行记录
                detectArray = obj.detect_array(channel, :);
                detectArrayNum = obj.detect_array_number(channel);
                detectArrayCount = obj.detect_array_count(channel, :);
                for pos_count = 1:length(binPos)
%                     disp("pos_count: " + pos_count);
                    arr = obj.detect_array(channel, 1:obj.detect_array_number(channel));
                    indices = obj.findBin(arr, binPos(pos_count), fft_x);
                    % indices = find(arr == binPos(pos_count) | abs(arr - binPos(pos_count)) == 1);
                    % is_member = ismember(obj.detect_array(channel, 1:obj.detect_array_number(channel)), binPos(pos_count));
                    % if obj.detect_array_number(channel) == 0 || sum(is_member) == 0
                    % 如果检测队列为空或者检测队列中没有该元素或没有与该元素绝对值相差1的值，直接加入
                    if obj.detect_array_number(channel) == 0 || isempty(indices)
                        % 检测队列的总数量加1
                        detectArrayNum = detectArrayNum + 1;
                        % 将峰对应的位置加入到检测队列中
                        detectArray(detectArrayNum) = binPos(pos_count);
                        % 并使其对应的count置1
                        detectArrayCount(detectArrayNum) = 1;
                    else % 如果检测队列有该元素或与该元素绝对值相差1的值，对其记录的count加1
                        % 找到对应的检测队列中的位置
                        % find_pos = find(obj.detect_array(channel, 1:obj.detect_array_number(channel)) == binPos(pos_count));
                        % 使其对应的count加1
                        detectArrayCount(indices) = detectArrayCount(indices) + 1;
                        % 修改纠正这个检测队列中的bin（调整1个bin）
                        detectArray(indices) = binPos(pos_count);
                    end
                end
                obj.detect_array(channel, :) = detectArray;
                obj.detect_array_number(channel) = detectArrayNum;
                obj.detect_array_count(channel, :) = detectArrayCount;

                % 将在检测队列中未能出现在此次preamble的峰值给剔除掉
                % 对检测队列中的每一个峰值进行遍历，找到没有在此次过滤出的峰出现的
                detectArrayNum = obj.detect_array_number(channel);
                countTmp = 1;
                for pos_count = 1:detectArrayNum
                    % 发现检测队里中的峰值在此次中未出现
                    binValue = obj.detect_array(channel, countTmp);
                    indices = obj.findBin(binPos, binValue, fft_x);
                    if isempty(indices)
                        % 在检测队列中进行删除
                        obj.detect_array(channel, countTmp:end) = [obj.detect_array(channel, countTmp+1:end), 0];
                        obj.detect_array_count(channel, countTmp:end) = [obj.detect_array_count(channel, countTmp+1:end), 0];
                        obj.detect_array_number(channel) = obj.detect_array_number(channel) - 1;
                    else
                        countTmp = countTmp + 1;
                    end
                end
                % 将检测队列中count数目等于7的元素放入检测成功的队列中
                detectArrayNum = obj.detect_array_number(channel);
                countTmp = 1;
                countThreshold = 5;
                for pos_count = 1:detectArrayNum
%                     if obj.detect_array_count(channel, countTmp) == 7
                    if obj.detect_array_count(channel, countTmp) == countThreshold
                        % 将挑选出的lora包的preamble bin记录
                        preambleBin = [preambleBin, obj.detect_array(channel, countTmp)];
                        % 将detect队列中对应元素删除，和上面做法相同
                        obj.detect_array(channel, countTmp:end) = [obj.detect_array(channel, countTmp+1:end), 0];
                        obj.detect_array_count(channel, countTmp:end) = [obj.detect_array_count(channel, countTmp+1:end), 0];
                        obj.detect_array_number(channel) = obj.detect_array_number(channel) - 1;
                    else
                        countTmp = countTmp + 1;
                    end
                end
            end
            % TODO: 如果没有active则清除队列
        end

        %{
            检测该信道该窗口内是否存在信号
            signals：已拆分信道和规定窗口后的信号
            isActive：布尔值，代表是否存在信号
        %}
        function [isActive] = isActive(obj, signals)
            dine = obj.loraSet.dine;
            fft_x = obj.loraSet.fft_x;
            isActive = false;
            % 通过FFT的峰值判断是否存在信号
            dechirp = signals .* obj.idealDownchirp;
            dechirp_fft = abs(fft(dechirp, dine));
            dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine);
            % 找到最大值
            [amp, ~] = max(dechirp_fft);
            % 最大值超过窗口均值的十倍
            if amp > mean(dechirp_fft) * 10
                isActive = true;
            end
        end

        %{
            根据loraSet中关于信道的信息，对输入的signals进行信道划分
            signals：需要拆分的原始信号
            splitSignal：拆分信道后的信号矩阵
        %}
        function [splitSignal] = divideChannel(obj, signals)
            splitSignal = [];
            for channel = 1 : obj.channelNum
                % TODO：将函数封装到类中
                signalOut = obj.signalFrequencyShift(signals, -obj.channelList(channel));
                signalOut = obj.lowPassFilterFir(signalOut);
                splitSignal = [splitSignal; signalOut];
            end
        end

        function binIndex = findBin(obj, binPos, findBin, fftX)
            arrLength = length(binPos);
            binPosRing = [binPos - fftX, binPos, binPos + fftX];
            index = find(binPosRing == findBin | abs(binPosRing - findBin) == 1, 1);
            % [~, closestIndex] = min(abs(binPosRing - findBin));
            if index <= arrLength
                binIndex = index;
            elseif index > arrLength*2
                binIndex = index - 2*arrLength;
            else
                binIndex = index - arrLength;
            end
        end
    end
end