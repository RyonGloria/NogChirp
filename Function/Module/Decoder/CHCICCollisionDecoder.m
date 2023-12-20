classdef CHCICCollisionDecoder < LoraDecoder
    properties
        signals;
        record;
        payloadBin;
        splitSignal;
        channelList;
        subchirpNum;
        preambleChannel;
        preambleSignal;
        preambleStartPos;
        preambleNum;
        preambleBin;
        channelArray;
        timeOffset;
        preambleEndPos;
        SFDPos;
        channelNum;
        channelMatrix;
        OverBandBw;
        window;
        errorFlag;
        errorMsg;
        peakStandard;
    end

    methods
        function obj = CHCICCollisionDecoder(loraSet, channelMatrix)
            obj@LoraDecoder(loraSet);
            obj.record = [];
            obj.payloadBin = [];
            
            % 获取配置中的一些参数
            bw = obj.loraSet.bw;
            OverBandBw = bw/2;     % 过渡带带宽
            obj.OverBandBw = OverBandBw;  % 记录过渡带带宽
            obj.channelNum = loraSet.channelNum;
            obj.subchirpNum = loraSet.subchirpNum;
            
            % 每个信道对应的中心频率（目前的信号已经进过频移至0频，其中心频率为采样时设置的中心频率）
            obj.channelList = zeros(1, obj.channelNum);
            for channel = 1:obj.channelNum
                % 记录每个信道的中心频率
                obj.channelList(channel) = (bw+OverBandBw)/2 - (obj.channelNum/2)*(bw+OverBandBw) ...
                    + (bw+OverBandBw)*(channel-1);
            end

            obj.channelMatrix = channelMatrix;
        end

        %{
            更改设置其中一些配置项
            注：对CHchirpCollisionDecoder重新初始化的时候建议使用该函数，避免重新创建该类
            channel：preamble所在信道
            window：对应的窗口位置
            preambleBin：preamble对应当前窗口的bin值
        %}
        function obj = setArgs(obj, channel, window, preambleBin)
            obj.preambleChannel = channel;
            obj.window = window;
            obj.preambleBin = preambleBin;
            obj.errorFlag = false; % 每次进行setArgs相当于一次初始化，所以要对errorFlag置为false
            obj.errorMsg = "";
            % 清空类中的某些中间值
            obj.splitSignal = [];
            obj.preambleSignal = [];
            obj.preambleEndPos = [];
            obj.cfo = [];
            obj.winOffset = [];
            obj.cfoDownchirp = [];
            obj.cfoUpchirp = [];
            obj.SFDPos = [];
            obj.channelArray = [];
            obj.timeOffset = [];
            obj.payloadBin = [];
        end

        function [obj] = decode(obj, signals)
            % 对信号进行信道划分
            obj = obj.divideChannel(signals);

            % 获取preamble所在信道的信号
            obj.preambleSignal = obj.splitSignal(obj.preambleChannel, :);

            % 找到对应信号的信号，检测preamble并找到最后一个preamble的位置
            obj = obj.detectPreambleEndPosBehind();
            if obj.errorFlag == true
                return;
            end
            
            % 通过preamble和SFD的bin来计算CFO和winoffset
            obj = obj.getcfoWinoff();
            if obj.errorFlag == true
                return;
            end

            % 调整信号的winoffset
            signals = circshift(signals, -round(obj.winOffset));
            obj.preambleSignal = circshift(obj.preambleSignal, -round(obj.winOffset));
           
            % 根据cfo重新生成带有decfo的idealchirp，用于解调
            obj = obj.rebuildIdealchirpCfo(0);

%             figure(3);
%             stft(signals(1:40*obj.loraSet.dine), obj.loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',obj.loraSet.fft_x);

            % 找到当前信号SFD的位置
            obj = obj.getSFDPos();
            if obj.errorFlag == true
                return;
            end
            
            % 获得payload阶段的跳信道矩阵(TODO: 不一定在窗口的14.25的位置处)
            obj = obj.getDownchirpSync();
            if obj.errorFlag == true
                return;
            end
            
            % 将信号四个信道滤波后保存在四个流里面
            obj = obj.divideChannel(signals);
            
            % 根据第一个跳信道subchirp（bin为0）来重新对齐
            obj = obj.alignWindows();
            if obj.errorFlag == true
                return;
            end
            
            % demodulate 解调信号
            % TODO：需要找出多个峰值，并采用投票的方式进行
            obj = obj.demodulate();
            if obj.errorFlag == true
                return;
            end
        end

        %{
            根据loraSet中关于信道的信息，对输入的signals进行信道划分
            signals：需要拆分的原始信号
            splitSignal：拆分信道后的信号矩阵
        %}
        function obj = divideChannel(obj, signals)
            obj.splitSignal = [];
            for channel = 1 : obj.channelNum
                signalOut = obj.signalFrequencyShift(signals, -obj.channelList(channel));
                signalOut = obj.lowPassFilterFir(signalOut);
                obj.splitSignal = [obj.splitSignal; signalOut];
            end
        end

        %{
            检测当前信道下对应premableBin的最后一个preamble的位置（为后续cfo计算准备）
            TODO: 因为涉及到多解码，建议把preambleEndPos信息作为输出传出
        %}
        function obj  = detectPreambleEndPosBehind(obj)
            dine = obj.loraSet.dine;
            fft_x = obj.loraSet.fft_x;
            preamble_len = obj.loraSet.Preamble_length;
            leakWidth = obj.loraSet.leakage_width1;
            window = obj.window;
            preambleSignal = obj.preambleSignal;
            preambleBin = obj.preambleBin;

            % 用于存储preambleLength+4数目窗口下，每个窗口的峰值
            candidate = cell(1, preamble_len + 4);
            % 处理找不到的情况
%             if window < 7
            if window < 4
                obj.errorFlag = true;
                obj.errorMsg = "detect发现的window窗口不满足要求";
                return;
            end
            % 搜寻第7个preamble前后的chirp
%             for t = window-7 : window+5
            for t = window-4 : window+8
                % 每一个窗口做FFT后，将获得的结果找出若干个峰值
                signal_tmp = preambleSignal(t*dine+1 : (t+1)*dine);
                dechirp = signal_tmp .* obj.idealDownchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine);
                % 找出若干个峰值，峰值间间隔为leakWidth
                result = obj.findpeaksWithShift(dechirp_fft, fft_x);
                % 记录峰值
%                 candidate{t - window + 8} = result(2, :);
                candidate{t - window + 5} = result(2, :);
            end

            % 找到每个窗口中，最接近preambleBin的值
            BinArray = obj.findClosetSyncWordBin(candidate, preambleBin);

            % 根据sync word的特征找到最后一个preamble
            Preamble_end_pos = obj.findPosWithSyncWord(BinArray);

            % 处理找不到的情况
            if exist('Preamble_end_pos', 'var') == 0 || Preamble_end_pos == 0
                obj.errorFlag = true;
                obj.errorMsg = "找不到最后一个preamble的位置";
                return;
            else
%                 obj.preambleEndPos = window - 7 + Preamble_end_pos - 2;
                obj.preambleEndPos = window - 4 + Preamble_end_pos - 2;
            end
        end

        %{
            从输入的元胞数组中，找到每一个元胞数组中与输入的bin接近的所有值
        %}
        function result = findClosetSyncWordBin(obj, cellArray, bin)
            fft_x = obj.loraSet.fft_x;
            % 找到每个元胞数组（每个元素都是一个数组）元素中，接近bin的所有值
            result = cell(1, length(cellArray));
            % 整理要查找的所有bin的可能值
            condition = zeros(3,3);
            tmp = -8;
            for i = 1:3 % 三种条件1,8,16
                tmp = tmp + 8;
                for k = 1:3 % 每种条件都为+1，-1，0的情况
                    condition(i, k) = bin + tmp + k - 2;
                    if condition(i, k) <= 0
                        condition(i, k) = condition(i, k) + fft_x;
                    elseif condition(i, k) > fft_x
                        condition(i, k) = condition(i, k) - fft_x;
                    end
                end
            end
            for i = 1:length(cellArray) % 遍历每一行
                row = cellArray{i}; % 获取当前行

                % 找到绝对值差值小于等于1的值的索引
                bin1 = row(ismember(row, condition(1, :)));
                % indices1 = find(abs(row - bin) <= 1);
                % 找到绝对值差值大于等于7小于等于9的值的索引
                bin2 = row(ismember(row, condition(2, :)));
                % indices2 = find(abs(row - bin) >= 7 & abs(row - bin) <= 9);
                % 找到绝对值差值大于等于15小于等于17的值的索引
                bin3 = row(ismember(row, condition(3, :)));
                % indices3 = find(abs(row - bin) >= 15 & abs(row - bin) <= 17);

                % 合并所有满足条件的索引
                allBin = unique([bin1, bin2, bin3]);

                % 如果找到了符合条件的数，存储它
                if ~isempty(allBin)
                    result{i} = allBin;
                else 
                    result{i} = 0;
                end
            end
        end

        %{
            在元胞数组（n*m的矩阵）中找到满足syncword特征的第一个窗口位置
            BinArray：元胞数组
            preamble_end_pos：满足sync word特征的
        %}
        function Preamble_end_pos = findPosWithSyncWord(obj, BinArray)
            fft_x = obj.loraSet.fft_x;
            Preamble_end_pos = 0;
            for loop1 = 3:length(BinArray)
                thirdArray = BinArray{loop1};
                secondArray = BinArray{loop1-1};
                firstArray = BinArray{loop1-2};
                for loop2 = 1:length(thirdArray)
                    thirdBin = thirdArray(loop2);
                    % 创建一个查找窗口
                    secondFindArr = zeros(1, 3);
                    for value = 7:9
                        secondFindArr(1, value-6) = thirdBin - value;
                        if secondFindArr(1, value-6) <= 0
                            secondFindArr(1, value-6) = secondFindArr(1, value-6) + fft_x;
                        end
                    end
                    for loop3 = 1:length(secondArray)
                        secondBin = secondArray(loop3);
                        if any(ismember(secondFindArr, secondBin))
                            firstFindArr = zeros(1, 3);
                            for value = 7:9
                                firstFindArr(1, value-6) = secondBin - value;
                                if firstFindArr(1, value-6) <= 0
                                    firstFindArr(1, value-6) = firstFindArr(1, value-6) + fft_x;
                                end
                            end
                            for loop4 = 1:length(firstArray)
                                firstBin = firstArray(loop4);
                                if any(ismember(firstFindArr, firstBin))
                                    Preamble_end_pos = loop1;
                                    return;
                                end
                            end
                        end
                    end
                end
            end
        end

        %{
            获取preamble bin和downchirp bin，计算得到cfo和winoffset
            signals：已划分信道的信号
            preambleBin：已确定的preamble的Bin值
        %}
        function obj = getcfoWinoff(obj)
            % 计算主峰的CFO(需要补零操作)
            % 对Preamble阶段的FFT峰值进行排序，得到前filter的峰
            zeropadding_size = obj.loraSet.factor;                   % 设置补零的数量，这里的decim表示，补上decim-1倍窗口的零，计算FFT时一共是decim倍的窗口（decim+1）
            d_sf = obj.loraSet.sf;
            d_bw = obj.loraSet.bw;
            dine = obj.loraSet.dine;
            fft_x = obj.loraSet.fft_x;
            leakWidth = obj.loraSet.leakage_width1;
            Preamble_num = obj.loraSet.Preamble_length;
            downchirp = obj.idealDownchirp;
            upchirp = obj.idealUpchirp;
            % 前面已经获得的最后一个preamble的窗口位置
            preambleEndPos = obj.preambleEndPos;
            preambleBin = obj.preambleBin;
            preambleSignal = obj.preambleSignal;
            
            dine_zeropadding = dine*zeropadding_size*2^(10-d_sf);
            fft_x_zeropadding = fft_x*zeropadding_size*2^(10-d_sf);

            % 获取最后一个preamble窗口的若干个峰值，找到最接近preambleBin的峰
            samples = preambleSignal((preambleEndPos-1)*dine+1 : (preambleEndPos)*dine);  % 倒数第二个preamble
            samples_fft = abs(fft(samples .* downchirp, dine_zeropadding, 2));
            samples_fft_merge = samples_fft(1:fft_x_zeropadding) + samples_fft(dine_zeropadding-fft_x_zeropadding+1:dine_zeropadding);
            % 找出若干个峰值，峰值间间隔为leakWidth
            result = obj.findpeaksWithShift(samples_fft_merge, fft_x_zeropadding);
            binPos = result(2, :);
            % 找到与bin最接近的bin的索引
%             closest_idx = obj.findClosetBin(binPos, preambleBin*zeropadding_size*2^(10-d_sf), fft_x_zeropadding);
            [~, closest_idx] = min(abs(binPos - preambleBin*zeropadding_size*2^(10-d_sf))); % 找出最接近bin值的元素的索引
            upchirpBin = binPos(closest_idx); % 获得最接近preambleBin的bin作为用来对齐的upchirpBin
            upchirpPeak = result(1, closest_idx); % 获取用于判断能量的标准值
            
            % 已知downchirp的位置，得到downchirp的bin（默认downchirp窗口内不存在downcrhip间的冲突）
            % TODO：可能需要考虑downcrhip发生冲突的问题
            SFD_samples = preambleSignal((preambleEndPos+3)*dine+1:(preambleEndPos+4)*dine);  % 第一个和第二个downchirp之间
            SFD_samples_fft = abs(fft(SFD_samples .* upchirp, dine_zeropadding));
            samples_fft_merge = SFD_samples_fft(1:fft_x_zeropadding) + SFD_samples_fft(dine_zeropadding-fft_x_zeropadding+1:dine_zeropadding);
            % 找出若干个峰值，峰值间间隔为leakWidth
            result = obj.findpeaksWithShift(samples_fft_merge, fft_x_zeropadding);
            if isempty(result)
                obj.errorFlag = true;
                obj.errorMsg = "找不到匹配的downchirp峰值";
                return;
            end
            [~, closest_idx] = min(abs(result(1, :) - upchirpPeak)); % 找出峰值能量最接近的元素的索引
            downchirpBin = result(2, closest_idx);
            
            % 计算CFO和窗口偏移量
            if upchirpBin + downchirpBin < fft_x_zeropadding*0.5
                cfo_bin = upchirpBin + downchirpBin - 2;
                obj.cfo = -cfo_bin/2/fft_x_zeropadding * d_bw;
                obj.winOffset = (downchirpBin - upchirpBin) / 2^(11-d_sf);
            elseif upchirpBin + downchirpBin > fft_x_zeropadding*1.5
                cfo_bin = upchirpBin + downchirpBin - fft_x_zeropadding*2 - 2;
                obj.cfo = -cfo_bin/2/fft_x_zeropadding * d_bw;
                obj.winOffset = (downchirpBin - upchirpBin) / 2^(11-d_sf);
            else
                cfo_bin = upchirpBin + downchirpBin - fft_x_zeropadding - 2;
                obj.cfo = -cfo_bin/2/fft_x_zeropadding * d_bw;
                obj.winOffset = (fft_x_zeropadding - (upchirpBin - downchirpBin)) / 2^(11-d_sf);
            end
        end

        function binIndex = findClosetBin(obj, binPos, findBin, fftX)
            arrLength = length(binPos);
            binPosRing = [binPos - fftX, binPos, binPos + fftX];
            [~, closestIndex] = min(abs(binPosRing - findBin));
            if closestIndex <= arrLength
                binIndex = closestIndex;
            elseif closestIndex > arrLength*2
                binIndex = closestIndex - 2*arrLength;
            else
                binIndex = closestIndex - arrLength;
            end
        end

        %{
            因为信号窗口被调整了，所以需要重新寻找SFD的位置，保证正确率
            preambleSignal：调整窗口后已划分信道的信号
        %}
        function obj = getSFDPos(obj)
            fftX = obj.loraSet.fft_x;
            dine = obj.loraSet.dine;
            leakWidth = obj.loraSet.leakage_width1;
%             preamble_len = obj.loraSet.Preamble_length;
            preambleBin = obj.preambleBin;
            preambleEndPos = obj.preambleEndPos;
            preambleSignal = obj.preambleSignal;
            if preambleEndPos < 2
                obj.errorFlag = true;
                obj.errorMsg = "preambleEndPos位置无效";
                return;
            end

            % 用于存储preambleLength+4数目窗口下，每个窗口的峰值
            candidate = cell(1, 5);
            
%             FFT_plot(preambleSignal, obj.loraSet, obj.cfoDownchirp, 5);
            % 同样思路，找到这几个窗口内的前20峰值，然后找到与bin1最接近的峰，判断其规律是否与sync word相同
            for windows = 1 : 5
%                 signal = preambleSignal((windows-1+preambleEndPos-preamble_len)*dine+1 : (windows+preambleEndPos-preamble_len)*dine);
                signal = preambleSignal((windows-1+preambleEndPos-2)*dine+1 : (windows+preambleEndPos-2)*dine);
                dechirp = signal .* obj.cfoDownchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
                result = obj.findpeaksWithShift(dechirp_fft, fftX);

                % 记录峰值
                candidate{windows} = result(2, :);
            end

            % 找到每个窗口中，最接近preambleBin的值
            BinArray = obj.findClosetSyncWordBin(candidate, 1);

            % 根据sync word的特征找到最后一个preamble
            SFDPos = obj.findPosWithSyncWord(BinArray);
            
            % 处理找不到最后一个preamble，报错
            if exist('SFDPos', 'var') == 0 || SFDPos == 0
                obj.errorFlag = true;
                obj.errorMsg = "找不到SFD的位置";
                return;
            else 
                obj.SFDPos = SFDPos + preambleEndPos - 2;
                % 获取SFD的能量，为后面找downchirp sync做准备
                peak = zeros(1, 2);
                for t = 0:1
                    signal = preambleSignal((obj.SFDPos+t)*dine+1:(obj.SFDPos+t+1)*dine);
                    dechirp = signal .* obj.cfoUpchirp;
                    dechirp_fft = abs(fft(dechirp, dine));
                    dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
                    % 冲突情况下，需要找到downchirp bin为1的峰的峰值
                    result = obj.findpeaksWithShift(dechirp_fft, fftX);
                    binPos = result(2, :);
                    closest_idx = obj.findClosetBin(binPos, 1, fftX);
                    if isempty(closest_idx)
                        obj.errorFlag = true;
                        obj.errorMsg = "找不到SFD中满足要求的峰";
                        return;
                    end
                    % [~, closest_idx] = min(abs(binPos - 1)); % 找出bin最接近1的元素的索引
                    peak(t+1) = result(1, closest_idx);
                end
                obj.peakStandard = mean(peak);
            end
        end

        %{
            根据SFD的位置，找到后面带有信道矩阵信息的downchirp，解调出信道矩阵信息
            preambleSignal：调整窗口后已划分信道的信号
        %}
        function obj = getDownchirpSync(obj)
            fftX = obj.loraSet.fft_x;
            dine = obj.loraSet.dine;
            preambleSignal = obj.preambleSignal;

            % dechirp做FFT，获得最大峰值即为信道矩阵的信息
            signal = preambleSignal((obj.SFDPos+4.25)*dine+1:(obj.SFDPos+5.25)*dine);
            dechirp = signal .* obj.cfoUpchirp;
            dechirp_fft = abs(fft(dechirp, dine));
            dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
            % 先获取若干个峰值，再从若干个峰值中找到能量最接近峰
            result = obj.findpeaksWithShift(dechirp_fft, fftX);
            % 找到差值最小的元素的索引
            [~, index] = min(abs(result(1, :) - obj.peakStandard));
            downchirp_bin = result(2, index);
            if isempty(downchirp_bin)
                obj.errorFlag = true;
                obj.errorMsg = "找不到带信号矩阵的bin值";
                return;
            end

            % 根据downchirp bin从生成的随机矩阵中获取跳信道矩阵
            obj.channelArray = obj.channelMatrix(downchirp_bin, :);
        end

        %{
            精对齐，根据信号规定第一个chirp所有subchirp的bin都为1的特性来进行对齐
        %}
        function obj = alignWindows(obj)
            fftX = obj.loraSet.fft_x * obj.loraSet.factor;
            dine = obj.loraSet.dine;
            subchirpNum = obj.subchirpNum;
            % 找到第一个chirp窗口的位置，截取出信号
            splitSignal = obj.splitSignal(:, (obj.SFDPos+5.25)*dine+1:(obj.SFDPos+6.25)*dine);
            
            peakCandidate = cell(1, subchirpNum);
            binCandidate = cell(1, subchirpNum);
            % 对第一个chirp的所有subchirp解码，确定精对齐的数目
            for subchirpCount = 0:subchirpNum-1
                chirpIntegrated = splitSignal(obj.channelArray(subchirpCount+1), subchirpCount*dine/subchirpNum+1 : (subchirpCount+1)*dine/subchirpNum);
                dechirp = obj.cfoDownchirp(subchirpCount*dine/subchirpNum+1 : (subchirpCount+1)*dine/subchirpNum) .* chirpIntegrated;
                % 进行补零的FFT
                dechirp_fft = abs(fft(dechirp, dine*obj.loraSet.factor));
                dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(end-fftX+1:end);
                % 找出若干个峰值，峰值间间隔为leakWidth
                result = obj.findpeaksWithShift(dechirp_fft, fftX);
                peakCandidate{subchirpCount + 1} = result(1, :);
                binCandidate{subchirpCount + 1} = result(2, :);
            end
            % 获得投票的bin值结果
            binVote = obj.getVoteResultWith2(binCandidate, obj.loraSet.factor);
            binVote = ceil(binVote / obj.loraSet.factor);
            % 在最后一个subchirp的bin序列中找到与结果接近的bin并根据峰值能量挑选正确的bin
            % 记录最接近的索引值和峰值
            closetIndex = 0;
            closetPeakDiff = realmax;
            % 遍历 binVote 中的元素
            binCandidateEnd = ceil(binCandidate{subchirpNum} / obj.loraSet.factor);
            for i = 1:length(binVote)
                value = binVote(i);
                % 找到在 binCandidate 中与 value 绝对值相差1的元素的索引
                % index = find(abs(binCandidateEnd - value) <= 1);
                index = obj.findClosetBin(binCandidateEnd, value, fftX);
                % 如果找到匹配的索引，则将其与最接近的峰值做比较
                if ~isempty(index)
                    if abs(peakCandidate{subchirpNum}(index) - obj.peakStandard/subchirpNum) < closetPeakDiff
                        closetIndex = index;
                        closetPeakDiff = abs(peakCandidate{subchirpNum}(index) - obj.peakStandard/subchirpNum);
                    end
                end
            end
            if closetIndex == 0
                obj.errorFlag = true;
                obj.errorMsg = "精对齐阶段找不到对齐的bin";
                return;
            end
            % 目前经验之谈是最后一个subchirp的bin能够实现精对齐
            obj.timeOffset = binCandidate{subchirpNum}(closetIndex);
        end

        %{
            从一个1*n的元胞数组中（每个元素都是一个不定长的向量），找到每个元胞数组元素中都存在相同的元素
            matrix：需要进行统计的元胞数组
            result：需要的结果
        %}
        function result = getVoteResultWith2(obj, cellArray, divisor)
            fft_x = obj.loraSet.fft_x;
            % 将元胞数组中的每个数组合并成一个单一的 cell 数组
            allValues = cat(2, cellArray{:});
            % 使用 unique 函数获取唯一的值
            uniqueValues = unique(allValues, 'rows', 'stable');

            % 初始化一个存储满足条件的数的数组
            result = [];

            % 遍历唯一值
            for value = uniqueValues
                valueDiv = ceil(value / divisor);
                % 检查每一行是否包含当前值或与当前值相差1的值
                if valueDiv == 1  % 避免1和fft_x的极端情况
                    valueSub1 = fft_x;
                    valueSub2 = fft_x-1;
                elseif valueDiv == 2
                    valueSub1 = 1;
                    valueSub2 = fft_x;
                else
                    valueSub1 = valueDiv - 1;
                    valueSub2 = valueDiv - 2;
                end

                if valueDiv == fft_x
                    valueAdd1 = 1;
                    valueAdd2 = 2;
                elseif valueDiv == fft_x - 1
                    valueAdd1 = fft_x;
                    valueAdd2 = 1;
                else
                    valueAdd1 = valueDiv + 1;
                    valueAdd2 = valueDiv + 2;
                end
                count = 0;
                for rowIndex = 1:length(cellArray)
                    row = ceil(cellArray{rowIndex} / divisor);
                    row_contains_value = any(ismember(row, [valueDiv, valueAdd1, valueAdd2, valueSub1, valueSub2]));
                    if row_contains_value
                        count = count + 1;
                    end
                end
                
                % 如果每一行都包含当前值或与当前值相差1的值，则将其存储
                if count == obj.subchirpNum
                    result = [result, value];
                end
            end
        end

        %{
            投票机制的解subchirp bin
        %}
        function obj = demodulate(obj)
            chirpNum = obj.loraSet.payloadNum;
            subchirpNum = obj.subchirpNum;
            factor = obj.loraSet.factor;
            fftX = obj.loraSet.fft_x;
            dine = obj.loraSet.dine;
            % 记录最后的bin值
            chirpBin = zeros(chirpNum, 1);
            % 定位到第一个chirp的窗口
            offset = 2*factor;  % 由于需要搜索窗口，需要给信号前面补充一些信号元素
            splitSignal = obj.splitSignal(:, (obj.SFDPos+5.25)*dine-obj.timeOffset+1 - offset:end);
            for chirp_count = 0:chirpNum-1
                % 1.对该chirp下所有subchirp做dechirpFFT，获得每个窗口的一系列峰值
                peakCandidate = cell(1, subchirpNum);
                binCandidate = cell(1, subchirpNum);
                for subchirpCount = 0:subchirpNum-1
                    chirpIntegrated = splitSignal(obj.channelArray(chirp_count*subchirpNum+subchirpCount+1), ...
                        chirp_count*dine + subchirpCount*dine/subchirpNum+1 + offset : ... 
                        chirp_count*dine + (subchirpCount+1)*dine/subchirpNum + offset);
                    dechirp = obj.cfoDownchirp(subchirpCount*dine/subchirpNum+1 : (subchirpCount+1)*dine/subchirpNum) .* chirpIntegrated;
                    dechirp_fft = abs(fft(dechirp, dine));
                    dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
                    % 获取这个窗口内最大的bin
                    % TODO：多冲突版本，投票机制
                    % 找出若干个峰值，峰值间间隔为leakWidth
                    result = obj.findpeaksWithShift(dechirp_fft, fftX);
                    peakCandidate{subchirpCount + 1} = result(1, :);
                    binCandidate{subchirpCount + 1} = result(2, :);
                end
                % 2.对subchirpNum窗口的值进行投票，得到有可能为这个chirp的bin值
                binVote = obj.getVoteResult(binCandidate, 1);
                % 记录最接近的索引值和峰值
                closetBin = 0;
                closetPeakDiff = realmax;
                % 遍历 binVote 中的元素
                for i = 1:length(binVote)
                    value = binVote(i);
                    energyAll = 0;  % 记录能量总和
                    % 计算对应的bin会在哪个信道发生了跳变
                    binPart = fftX/subchirpNum;
                    jumpPart = subchirpNum - floor(value/binPart);
                    % 遍历binCandidate中的每一行
                    for t = 1:length(binCandidate)
                        if t == jumpPart  % 如果这个subchirp位置处发生了跳变，则不用来做能量判断
                            continue;
                        end
                        % 找到当前行与目标值绝对值相差1的元素的索引
                        % index = find(abs(binCandidate{t} - value) <= 1);
                        index = obj.findClosetBin(binCandidate{t}, value, fftX);
                        energyAll = energyAll + peakCandidate{t}(index);
                    end 
                    if abs(energyAll - obj.peakStandard*(1-1/subchirpNum)) < closetPeakDiff
                        closetBin = value;
                        closetPeakDiff = abs(energyAll - obj.peakStandard*(1-1/subchirpNum));
                    end
                end
                % 获得最有可能得bin值
                bin = closetBin;
                % 3.通过滑动窗口和取众数来确定这个正确的bin值
                binTmp = zeros(subchirpNum, factor+1); % 记录这个滑动窗口中与bin最接近的bin
                % 对每个chirp窗口进行滑动获取bin值
                for slideWindow = -factor/2:factor/2
                    for subchirpCount = 0:subchirpNum-1
                        chirpIntegrated = splitSignal(obj.channelArray(chirp_count*subchirpNum+subchirpCount+1), ...
                            chirp_count*dine + subchirpCount*dine/subchirpNum+1 + offset + slideWindow : ... 
                            chirp_count*dine + (subchirpCount+1)*dine/subchirpNum + offset + slideWindow);
                        dechirp = obj.cfoDownchirp(subchirpCount*dine/subchirpNum+1 : (subchirpCount+1)*dine/subchirpNum) .* chirpIntegrated;
                        dechirp_fft = abs(fft(dechirp, dine));
                        dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
                        % 找出若干个峰值，峰值间间隔为leakWidth
                        result = obj.findpeaksWithShift(dechirp_fft, fftX);
                        % index = find(abs(result(2, :) - bin) <= 1);
                        index = obj.findClosetBin(result(2, :), bin, fftX);
                        if isempty(result) || isempty(index)
                            obj.errorFlag = true;
                            obj.errorMsg = "解subchirpbin投票之后找不到峰值";
                            return;
                        end
                        binTmp(subchirpCount+1, slideWindow + factor/2 + 1) = result(2, index);
                    end
                end
                % 投票选出峰值
                chirpBin(chirp_count + 1) = mode(binTmp(:));
            end
            obj.payloadBin = chirpBin;
        end

        %{
            从一个1*n的元胞数组中（每个元素都是一个不定长的向量），找到每个元胞数组元素中都存在相同的元素
            matrix：需要进行统计的元胞数组
            result：需要的结果
        %}
        function result = getVoteResult(obj, cellArray, divisor)
            fft_x = obj.loraSet.fft_x;
            % 将元胞数组中的每个数组合并成一个单一的 cell 数组
            allValues = cat(2, cellArray{:});
            % 使用 unique 函数获取唯一的值
            uniqueValues = unique(allValues, 'rows', 'stable');

            % 初始化一个存储满足条件的数的数组
            result = [];

            % 遍历唯一值
            for value = uniqueValues
                valueDiv = ceil(value / divisor);
                % 检查每一行是否包含当前值或与当前值相差1的值
                if valueDiv == 1  % 避免1和fft_x的极端情况
                    valueSub1 = fft_x;
                else
                    valueSub1 = valueDiv - 1;
                end

                if valueDiv == fft_x
                    valueAdd1 = 1;
                else
                    valueAdd1 = valueDiv + 1;
                end
                count = 0;
                for rowIndex = 1:length(cellArray)
                    row = ceil(cellArray{rowIndex} / divisor);
                    row_contains_value = any(ismember(row, [valueDiv, valueAdd1, valueSub1]));
                    if row_contains_value
                        count = count + 1;
                    end
                end
                
                % 如果每一行都包含当前值或与当前值相差1的值，则将其存储
                if count == obj.subchirpNum
                    result = [result, value];
                end
            end
        end

    end
end