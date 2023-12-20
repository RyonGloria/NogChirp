classdef CHchirpCollisionDecoder < LoraDecoder
    properties
        payloadBin;
        splitSignal;
        channelList;
        subchirpNum;
        preambleChannel;
        preambleSignal;
        preambleStartPos;
        Downchirp_ind;
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
        SFDPeakAmp;
        DebugUtil;
        preamblePeak;
    end

    methods
        function obj = CHchirpCollisionDecoder(loraSet, channelMatrix, DebugUtil)
            obj@LoraDecoder(loraSet);
            % 初始化debug工具
            obj.DebugUtil = DebugUtil;
            obj.DebugUtil.info("", "初始化CHchirpCollisionDecoder");

            % 获取配置中的一些参数
            bw = obj.loraSet.bw;
            OverBandBw = bw/2;     % 过渡带带宽
            obj.OverBandBw = OverBandBw;  % 记录过渡带带宽
            obj.DebugUtil.debug("\t", "过渡带：" + obj.OverBandBw);
            obj.channelNum = loraSet.channelNum;
            obj.subchirpNum = loraSet.subchirpNum;
            obj.DebugUtil.debug("\t", "信道数目：" + obj.channelNum);
            obj.DebugUtil.debug("\t", "subchirp数目：" + obj.subchirpNum);
            
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
            obj.DebugUtil.info("", "根据detect获取的参数，重新初始化decoder的参数");
            obj.preambleChannel = channel;
            obj.DebugUtil.debug("\t", "检测到preamble信道位置：" + obj.preambleChannel);
            obj.window = window;
            obj.DebugUtil.debug("\t", "检测到preamble窗口位置：" + obj.window);
            obj.preambleBin = preambleBin;
            obj.DebugUtil.debug("\t", "检测到preamble对应的bin：" + obj.preambleBin);
        end

        function obj = clear(obj)
            % 清空类中的某些中间值
            obj.DebugUtil.info("", "清空decode中间过程所有属性值");
            obj.payloadBin = [];
            obj.splitSignal = [];
            obj.preambleSignal = [];
            obj.preambleStartPos = [];
            obj.channelArray = [];
            obj.timeOffset = [];
            obj.preambleEndPos = [];
            obj.SFDPos = [];
            obj.peakStandard = [];
            obj.SFDPeakAmp = [];
            obj.cfo = [];
            obj.winOffset = [];
            obj.cfoDownchirp = [];
            obj.cfoUpchirp = [];
            obj.preamblePeak = [];
            obj.errorFlag = false; % 每次进行setArgs相当于一次初始化，所以要对errorFlag置为false
            obj.errorMsg = "";
        end

        function [obj] = decode(obj, signals)
            obj.DebugUtil.info("", "调用CHchirpCollisionDecoder=====================================");
            % 清除所有中间变量属性
            obj = obj.clear();

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
            
%             obj = obj.falsePositiveJudge();
%             if obj.errorFlag == true
%                 return;
%             end
            
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
            obj.DebugUtil.info("", "CHchirpCollisionDecoder运行结束======================================");
        end

        function obj = detectPass(obj, signals)
            obj.DebugUtil.info("", "调用CHchirpCollisionDecoder,detectPass=====================================");
            % 清除所有中间变量属性
            obj = obj.clear();

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

            % 找到当前信号SFD的位置
            obj = obj.getSFDPos();
            if obj.errorFlag == true
                return;
            end
            obj.DebugUtil.info("", "CHchirpCollisionDecoder,detectPass运行结束======================================");
        end

        function [obj] = decodeDebug(obj, signals)
            obj.DebugUtil.info("", "调用CHchirpCollisionDecoder=====================================");
            % 清除所有中间变量属性
            obj = obj.clear();

            % 对信号进行信道划分
            obj = obj.divideChannel(signals);

            % 获取preamble所在信道的信号
            obj.preambleSignal = obj.splitSignal(obj.preambleChannel, :);

            obj = obj.corrDownchirp(obj.window, 10);
            if obj.errorFlag == true
                return;
            end

            obj = obj.getPrembleEndPos();

            % figure(1);
            % FFT_plot(obj.preambleSignal((obj.preambleEndPos-1)*obj.loraSet.dine+1:end), obj.loraSet, obj.idealDownchirp, 6);

            obj = obj.getcfoWinoff();
            if obj.errorFlag == true
                return;
            end

            % 调整信号的winoffset
            signals = circshift(signals, -round(obj.winOffset));
            obj.preambleSignal = circshift(obj.preambleSignal, -round(obj.winOffset));
           
            % 根据cfo重新生成带有decfo的idealchirp，用于解调
            obj = obj.rebuildIdealchirpCfo(0);

            % 找到当前信号SFD的位置
            % obj = obj.getSFDPos();
            obj = obj.corrDownchirp(obj.preambleEndPos - 2, 10);
            if obj.errorFlag == true
                return;
            end

            obj = obj.getSFDPosCIC();
%             figure(1);
%             dine = obj.loraSet.dine;
%             FFT_plot(signals((obj.SFDPos-1)*dine+1:end), obj.loraSet, obj.cfoUpchirp, 6);
            
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
            obj.DebugUtil.info("", "CHchirpCollisionDecoder运行结束======================================");
        end

        function [obj] = corrDownchirp(obj, window, windowNum)
            obj.Downchirp_ind = [];
            % window = obj.window;
            dine = obj.loraSet.dine;
            fft_x = obj.loraSet.fft_x;
            pnts_threshold = 40;      % Max. # of peaks to extract from Corrrelation Plot
            num_DC = 2.25;
            preambleSignal = obj.preambleSignal(window*dine+1 : (window+windowNum)*dine);
            preambleSignal = downsample(preambleSignal, obj.loraSet.factor);
            downchirp = obj.idealDownchirp;
            downchirp = downsample(downchirp, obj.loraSet.factor);
            % Cross Correlation with a Single downchirp
            Cross_Corr = zeros(1, length(preambleSignal) - fft_x - 1);
            for i = 1:length(preambleSignal) - fft_x - 1
                Cross_Corr(i) = sum(preambleSignal( i : i + fft_x - 1) .* conj(downchirp))...
                        / sqrt(sum( preambleSignal( i : i + fft_x - 1) .* conj(preambleSignal( i : i + fft_x - 1)) ) * ...
                        sum( downchirp .* conj(downchirp)));
            end
            Cross_Corr = Cross_Corr(isfinite(Cross_Corr));
            corr_threshold =  4*sum(abs(Cross_Corr))/length(Cross_Corr);
%             plot(abs(Cross_Corr));

            n_samp_array = [];
            peak_ind_prev = [];
            for i = 0:floor(length(Cross_Corr)/fft_x)-1
                % windowing Cross-Correlation (window length fft_x samples)
                wind = abs(Cross_Corr(i*fft_x + 1 : (i+1) * fft_x));                            
                % Extract Multiple Correlation Peaks
                peak_ind_curr = obj.get_max(wind,corr_threshold,pnts_threshold);          
                if(length(peak_ind_prev) ~= 0 && length(peak_ind_curr) ~= 0)
                    for j = 1:length(peak_ind_curr)
                        for k = 1:length(peak_ind_prev)
                            % check if combination of any two peaks in consecutive window are N samples apart
                            if(peak_ind_curr(j) == peak_ind_prev(k))                    
                                n_samp_array = [n_samp_array  peak_ind_prev(k)+((i-1)*fft_x) peak_ind_curr(j)+((i)*fft_x)];
                            end
                            % This extracts a list of all peaks that are N samples
                            % apart
                        end
                    end 
                end
                peak_ind_prev = peak_ind_curr;
            end

            Downchirp_ind = [];
            for i = 1:length(n_samp_array)
                c = 0;
                ind_arr = n_samp_array(i) : fft_x : n_samp_array(i) + (fft_x);
                
                for j = 1:length(ind_arr)
                    c = c + sum( n_samp_array == ind_arr(j) );
                end
                % Find from the list all the peaks that appear consecutively for
                % more than 2 windows (Downchirp should give 3 peaks, N sampled apart)
                if( c >= 2 )
                    Downchirp_ind = [Downchirp_ind; [ind_arr]];
                end
            end

            % filter Downchirps that are with in 3 samples (same pkt detected twice due to peak energy spread)
            temp = [];
            indices = [zeros(1,floor(num_DC)); Downchirp_ind];
            for i = 2:size(indices,1)
                
                if(isempty(temp))
                    temp = [temp; indices(i,:)];
                else
                    if( min(abs(indices(i) - temp(:,1))) > 3 )
                        temp = [temp; indices(i,:)];
                    end
                end
            end
            Downchirp_ind = temp;
            if isempty(Downchirp_ind)
                obj.DebugUtil.warning("\t\t", "对齐前找不到SFD的位置");
                obj.errorFlag = true;
                obj.errorMsg = "对齐前找不到SFD的位置";
                return;
            else
                obj.Downchirp_ind = Downchirp_ind;
            end
        end

        function obj = getPrembleEndPos(obj)
            fft_x = obj.loraSet.fft_x;
            obj.preambleEndPos = obj.window + round(obj.Downchirp_ind(1,1) / fft_x) - 2;
        end

        function obj = getSFDPosCIC(obj)
            dine = obj.loraSet.dine;
            fft_x = obj.loraSet.fft_x;
            obj.SFDPos = obj.preambleEndPos - 2 + round(obj.Downchirp_ind(1,1) / fft_x);
            obj.SFDPeakAmp = zeros(1, 2);
            for t = 0:1
                signal = obj.preambleSignal((obj.SFDPos+t)*dine+1:(obj.SFDPos+t+1)*dine);
                dechirp = signal .* obj.cfoUpchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine);
                % 冲突情况下，需要找到downchirp bin为1的峰的峰值
                result = obj.findpeaksWithShift(dechirp_fft, fft_x);
                binPos = result(2, :);
                closest_idx = obj.findClosetBin(binPos, 1, fft_x);
                if isempty(closest_idx)
                    obj.errorFlag = true;
                    obj.errorMsg = "找不到SFD中满足要求的峰";
                    return;
                end
                peak(t+1) = result(1, closest_idx);
                obj.SFDPeakAmp(t+1) = peak(t+1);
            end
            obj.peakStandard = mean(peak);
        end

        function [out] = get_max(obj, arr, threshold, num_pnts)
        %GET_MAX, this function returns uptil num_pnts many maximums above
            out = [];
            for i = 1:num_pnts
                [a, b] = max(arr);
                if(a < threshold)
                    return;
                else
                    out = [out b];
                    arr(b) = 0;
                end
            end
        
        end

        function obj = falsePositiveJudge(obj)
            fftX = obj.loraSet.fft_x;
            dine = obj.loraSet.dine;
            preambleSignal = obj.preambleSignal;

            % preamble peak
            preamblePeak = 0;
            for index = 0:1
                signal = preambleSignal((obj.SFDPos-4+index)*dine+1:(obj.SFDPos-3+index)*dine);
                dechirp = signal .* obj.cfoDownchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
                preamblePeak = preamblePeak + dechirp_fft(1);
            end
            preamblePeak = preamblePeak/2;

            % sync word peak
            syncWordPeak = 0;
            for index = 0:1
                signal = preambleSignal((obj.SFDPos-2+index)*dine+1:(obj.SFDPos-1+index)*dine);
                dechirp = signal .* obj.cfoDownchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
                syncWordPeak = syncWordPeak + dechirp_fft(index*8+9);
            end
            syncWordPeak = syncWordPeak/2;


            % SFD peak
            sfdPeak = 0;
            for index = 0:1
                signal = preambleSignal((obj.SFDPos+index)*dine+1:(obj.SFDPos+1+index)*dine);
                dechirp = signal .* obj.cfoUpchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
                sfdPeak = sfdPeak + dechirp_fft(1);
            end
            sfdPeak = sfdPeak/2;

            thresholdMin = 0.9;
            thresholdMax = 1.1;
            condition1 = preamblePeak/syncWordPeak > thresholdMin && preamblePeak/syncWordPeak < thresholdMax;
            condition2 = syncWordPeak/sfdPeak > thresholdMin && syncWordPeak/sfdPeak < thresholdMax;
            condition3 = sfdPeak/preamblePeak > thresholdMin && sfdPeak/preamblePeak < thresholdMax;
            if ~(condition1 && condition2 && condition3)
                obj.errorFlag = true;
                obj.errorMsg = "检测为假阳性（非正常LoRa包）";
                return;
            end
        end

        %{
            根据loraSet中关于信道的信息，对输入的signals进行信道划分
            signals：需要拆分的原始信号
            splitSignal：拆分信道后的信号矩阵
        %}
        function obj = divideChannel(obj, signals)
            obj.DebugUtil.info("\t", "对信道进行拆分");
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
            obj.DebugUtil.info("\t", "找到最后一个preamble的窗口位置");
            dine = obj.loraSet.dine;
            fft_x = obj.loraSet.fft_x;
            preamble_len = obj.loraSet.Preamble_length;
            leakWidth = obj.loraSet.leakage_width1;
            window = obj.window;
            preambleSignal = obj.preambleSignal;
            preambleBin = obj.preambleBin;

            % 用于存储preambleLength+4数目窗口下，每个窗口的峰值
            candidate = cell(1, preamble_len);
            cancidateAmp = cell(1, preamble_len);
            % 处理找不到的情况
%             if window < 7
            if window > 50
                obj.DebugUtil.warning("\t\t", "preamble输入窗口不满足要求，窗口位置为" + window);
                obj.errorFlag = true;
                obj.errorMsg = "detect发现的window窗口不满足要求";
                return;
            end
            % 搜寻第7个preamble前后的chirp
            for t = window : window+7
                % 每一个窗口做FFT后，将获得的结果找出若干个峰值
                signal_tmp = preambleSignal(t*dine+1 : (t+1)*dine);
                dechirp = signal_tmp .* obj.idealDownchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine);
                % 找出若干个峰值，峰值间间隔为leakWidth
                result = obj.findpeaksWithShift(dechirp_fft, fft_x);
                % 记录峰值
                candidate{t - window + 1} = result(2, :);
                cancidateAmp{t - window + 1} = result(1, :);
            end

            % 找到每个窗口中，最接近preambleBin的值
            [BinArray, AmpArray] = obj.findClosetSyncWordBin(candidate, cancidateAmp, preambleBin);
            obj.DebugUtil.debug("\t\t", "找到窗口中，最接近preambleBin的值");
            obj.DebugUtil.debug("\t\t", "BinArray");
            obj.DebugUtil.debug("\t\t", BinArray);

            % 根据sync word的特征找到最后一个preamble
            Preamble_end_pos = obj.findPosWithSyncWord(BinArray, AmpArray);
            obj.DebugUtil.debug("\t\t", "根据sync word的特征找到最后一个preamble");
            obj.DebugUtil.debug("\t\t", "Preamble_end_pos：");
            obj.DebugUtil.debug("\t\t", Preamble_end_pos);

            % 处理找不到的情况
            if exist('Preamble_end_pos', 'var') == 0 || Preamble_end_pos == 0
                obj.DebugUtil.warning("\t\t", "找不到最后一个preamble的位置");
                obj.errorFlag = true;
                obj.errorMsg = "找不到最后一个preamble的位置";
                return;
            else
                obj.preambleEndPos = window + Preamble_end_pos - 2;
                obj.DebugUtil.info("\t", "找到最后一个preamble的窗口位置：" + obj.preambleEndPos);
                % 获取preamble的能量
                signal_tmp = preambleSignal((obj.preambleEndPos-1)*dine+1 : (obj.preambleEndPos)*dine);
                dechirp = signal_tmp .* obj.idealDownchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine);
                % 找出若干个峰值，峰值间间隔为leakWidth
                result = obj.findpeaksWithShift(dechirp_fft, fft_x);
                binPos = result(2, :);
                % 记录峰值
                [~, closest_idx] = min(abs(binPos - preambleBin));
                obj.preamblePeak = result(1, closest_idx);
            end
        end

        %{
            从输入的元胞数组中，找到每一个元胞数组中与输入的bin接近的所有值
        %}
        function [result, ampResult] = findClosetSyncWordBin(obj, cellArray, ampArray, bin)
            fft_x = obj.loraSet.fft_x;
            % 找到每个元胞数组（每个元素都是一个数组）元素中，接近bin的所有值
            result = cell(1, length(cellArray));
            ampResult = cell(1, length(cellArray));
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
                rowAmp = ampArray{i};

                % 找到绝对值差值小于等于1的值的索引
                bin1 = row(ismember(row, condition(1, :)));
                amp1 = rowAmp(ismember(row, condition(1, :)));
                % indices1 = find(abs(row - bin) <= 1);
                % 找到绝对值差值大于等于7小于等于9的值的索引
                bin2 = row(ismember(row, condition(2, :)));
                amp2 = rowAmp(ismember(row, condition(2, :)));
                % indices2 = find(abs(row - bin) >= 7 & abs(row - bin) <= 9);
                % 找到绝对值差值大于等于15小于等于17的值的索引
                bin3 = row(ismember(row, condition(3, :)));
                amp3 = rowAmp(ismember(row, condition(3, :)));
                % indices3 = find(abs(row - bin) >= 15 & abs(row - bin) <= 17);

                % 合并所有满足条件的索引
                allBin = [bin1, bin2, bin3];
                allAmp = [amp1, amp2, amp3];

                % 如果找到了符合条件的数，存储它
                if ~isempty(allBin)
                    result{i} = allBin;
                    ampResult{i} = allAmp;
                else 
                    result{i} = 0;
                    ampResult{i} = 0;
                end
            end
        end

        %{
            在元胞数组（n*m的矩阵）中找到满足syncword特征的第一个窗口位置
            BinArray：元胞数组
            preamble_end_pos：满足sync word特征的
        %}
        function Preamble_end_pos = findPosWithSyncWord(obj, BinArray, AmpArray)
            fft_x = obj.loraSet.fft_x;
            Preamble_end_pos = 0;
            for loop1 = 3:length(BinArray)
                thirdArray = BinArray{loop1};
                thirdAmpArray = AmpArray{loop1};
                secondArray = BinArray{loop1-1};
                secondAmpArray = AmpArray{loop1-1};
                firstArray = BinArray{loop1-2};
                firstAmpArray = AmpArray{loop1-2};
                for loop2 = 1:length(thirdArray)
                    thirdBin = thirdArray(loop2);
                    thirdAmp = thirdAmpArray(loop2);
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
                        secondAmp = secondAmpArray(loop3);
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
                                firstAmp = firstAmpArray(loop4);
                                if any(ismember(firstFindArr, firstBin))
                                    if thirdAmp/firstAmp >= 0.5 && secondAmp/thirdAmp >= 0.5 && secondAmp/thirdAmp <= 1.5
                                        Preamble_end_pos = loop1;
                                        return;
                                    end
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
            obj.DebugUtil.info("\t", "getcfoWinoff，获取cfo和winOff");
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
            % [~, closest_idx] = min(abs(binPos - preambleBin*zeropadding_size*2^(10-d_sf))); % 找出最接近bin值的元素的索引
            % 找到与bin范围内接近的所有峰
            binArr = obj.findBinRange(binPos, preambleBin*zeropadding_size*2^(10-d_sf), zeropadding_size*2, fft_x_zeropadding);
            [~, closest_idx] = min(abs(samples_fft_merge(binArr) - obj.preamblePeak)); 
            % upchirpBin = binPos(closest_idx); % 获得最接近preambleBin的bin作为用来对齐的upchirpBin
            upchirpBin = binArr(closest_idx);
            upchirpPeak = samples_fft_merge(upchirpBin); % 获取用于判断能量的标准值
            obj.DebugUtil.debug("\t\t", "upchirpBin:" + string(upchirpBin));
            obj.DebugUtil.debug("\t\t", "upchirpPeak:" + string(upchirpPeak));
            
            % 已知downchirp的位置，得到downchirp的bin（默认downchirp窗口内不存在downcrhip间的冲突）
            % TODO：可能需要考虑downcrhip发生冲突的问题
            SFD_samples = preambleSignal((preambleEndPos+3)*dine+1:(preambleEndPos+4)*dine);  % 第一个和第二个downchirp之间
            SFD_samples_fft = abs(fft(SFD_samples .* upchirp, dine_zeropadding));
            samples_fft_merge = SFD_samples_fft(1:fft_x_zeropadding) + SFD_samples_fft(dine_zeropadding-fft_x_zeropadding+1:dine_zeropadding);
            % 找出若干个峰值，峰值间间隔为leakWidth
            result = obj.findpeaksWithShift(samples_fft_merge, fft_x_zeropadding);
            if isempty(result)
                obj.DebugUtil.warning("", "找不到匹配的downchirp峰值");
                obj.errorFlag = true;
                obj.errorMsg = "找不到匹配的downchirp峰值";
                return;
            end
            [~, closest_idx] = min(abs(result(1, :) - upchirpPeak)); % 找出峰值能量最接近的元素的索引
            downchirpBin = result(2, closest_idx);
            obj.DebugUtil.debug("\t\t", "downchirpBin:" + string(downchirpBin));
            
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
            obj.DebugUtil.info("\t", "获得CFO：" + obj.cfo);
            obj.DebugUtil.info("\t", "获得winOffset：" + obj.winOffset);
        end

        function binPosRing = findBinRange(obj, binPos, findBin, range, fftX)
            binPosRing = [binPos - fftX, binPos, binPos + fftX];
            binPosRing = binPosRing(binPosRing >= findBin - range & binPosRing <= findBin + range);
            % [~, closestIndex] = min(abs(binPosRing - findBin));
            for index = 1:length(binPosRing)
                if binPosRing(index) <= 0
                    binPosRing(index) = binPosRing(index) + fftX;
                elseif binPosRing(index) > fftX
                    binPosRing(index) = binPosRing(index) - fftX;
                end
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
            obj.DebugUtil.info("\t", "getSFDPos，获取SFD的位置")
            fftX = obj.loraSet.fft_x;
            dine = obj.loraSet.dine;
            leakWidth = obj.loraSet.leakage_width1;
%             preamble_len = obj.loraSet.Preamble_length;
            preambleBin = obj.preambleBin;
            preambleEndPos = obj.preambleEndPos;
            preambleSignal = obj.preambleSignal;
            if preambleEndPos < 2
                obj.DebugUtil.warning("\t", "preambleEndPos位置无效");
                obj.errorFlag = true;
                obj.errorMsg = "preambleEndPos位置无效";
                return;
            end

            % 用于存储preambleLength+4数目窗口下，每个窗口的峰值
            candidate = cell(1, 5);
            candidatePeaks = cell(1, 5);
            
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
                candidatePeaks{windows} = result(1, :);
            end

            % 找到每个窗口中，最接近preambleBin的值
            [BinArray, AmpArray] = obj.findClosetSyncWordBin(candidate, candidatePeaks, 1);
            obj.DebugUtil.debug("\t\t", "找到每个窗口中，最接近preambleBin的值");
            obj.DebugUtil.debug("\t\t", "BinArray：");
            obj.DebugUtil.debug("\t\t", BinArray);

            % 根据sync word的特征找到最后一个preamble
            SFDPos = obj.findPosWithSyncWord(BinArray, AmpArray);
            
            
            % 处理找不到最后一个preamble，报错
            if exist('SFDPos', 'var') == 0 || SFDPos == 0
                obj.DebugUtil.warning("\t\t", "找不到SFD的位置");
                obj.errorFlag = true;
                obj.errorMsg = "找不到SFD的位置";
                return;
            else 
                obj.SFDPos = SFDPos + preambleEndPos - 2;
                obj.DebugUtil.info("\t\t", "根据sync word的特征找到第一个SFD");
                obj.DebugUtil.info("\t\t", "SFDPos: " + obj.SFDPos);
                % 获取SFD的能量，为后面找downchirp sync做准备
                peak = zeros(1, 2);
                obj.SFDPeakAmp = zeros(1, 2);
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
                    peak(t+1) = result(1, closest_idx);
                    obj.SFDPeakAmp(t+1) = peak(t+1);
                end
                obj.peakStandard = mean(peak);
                obj.DebugUtil.debug("\t\t", "SFD的两个峰值" + obj.SFDPeakAmp(1) + ", " + obj.SFDPeakAmp(2));
                obj.DebugUtil.debug("\t\t", "SFD的峰值均值" + obj.peakStandard);
            end
        end

        %{
            在元胞数组（n*m的矩阵）中找到满足syncword特征的第一个窗口位置
            BinArray：元胞数组
            preamble_end_pos：满足sync word特征的
        %}
        function Preamble_end_pos = findPosWithSyncWordWithAmp(obj, BinArray, preambleBin)
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
            根据SFD的位置，找到后面带有信道矩阵信息的downchirp，解调出信道矩阵信息
            preambleSignal：调整窗口后已划分信道的信号
        %}
        function obj = getDownchirpSync(obj)
            obj.DebugUtil.info("\t", "getDownchirpSync,获取跳信道矩阵信息");
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
            obj.DebugUtil.info("\t\t", "downchirp_bin: " + downchirp_bin);
            if isempty(downchirp_bin)
                obj.DebugUtil.warning("\t", "找不到带信号矩阵的bin值");
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
            obj.DebugUtil.info("\t", "alignWindows, 精对齐窗口");
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
            obj.DebugUtil.debug("\t\t", "投票所有可能结果：");
            obj.DebugUtil.debug("\t\t", "binVote");
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
                obj.DebugUtil.warning("\t", "精对齐阶段找不到对齐的bin");
                obj.errorFlag = true;
                obj.errorMsg = "精对齐阶段找不到对齐的bin";
                return;
            end
            % 目前经验之谈是最后一个subchirp的bin能够实现精对齐
            obj.timeOffset = binCandidate{subchirpNum}(closetIndex);
            obj.DebugUtil.info("\t", "得到精对齐的timeOffset：" + obj.timeOffset);
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
            obj.DebugUtil.info("\t", "demodulate，解调出payloadBin");
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
            obj.DebugUtil.info("\t", "获取所有的的paylodaBin: ");
            obj.DebugUtil.info("\t", chirpBin');
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