classdef CHchirpDecoder < LoraDecoder
    properties
        record;
        payloadBin;
        splitSignal;
        channelList;
        subchirpNum;
        preambleActiveChannel;
        preambleSignal;
        preambleStartPos;
        preambleNum;
        channelArray;
        timeOffset;
        preambleEndPos;
        SFDPos;
        channelNum;
        channelMatrix;
        PreambleBin;
        upchirpbin;
        downchirpbin;
    end
    
    methods
        % 初始化方法
        function obj = CHchirpDecoder(loraSet)
            obj@LoraDecoder(loraSet);
            obj.record = [];
            obj.payloadBin = [];
            
            % 获取配置中的一些参数
            bw = obj.loraSet.bw;
            OverBandBw = bw/2;     % 过渡带带宽
            obj.record.OverBandBw = OverBandBw;  % 记录过渡带带宽
            obj.channelNum = loraSet.channelNum;
            obj.subchirpNum = loraSet.subchirpNum;
            
            % 每个信道对应的中心频率（目前的信号已经进过频移至0频，其中心频率为采样时设置的中心频率）
            obj.channelList = zeros(1, obj.channelNum);
            for channel = 1:obj.channelNum
                % 记录每个信道的中心频率
                obj.channelList(channel) = (bw+OverBandBw)/2 - (obj.channelNum/2)*(bw+OverBandBw) ...
                    + (bw+OverBandBw)*(channel-1);
            end

            obj.channelMatrix = [];
            for i = 0:obj.loraSet.fft_x-1
                obj.channelMatrix = [obj.channelMatrix; obj.generateArray(i, obj.channelNum)];
            end
        end
        
        function obj = decode(obj, signals)
            
            % 将信号拆分成channelNum个信道信号
            obj = obj.divideChannel(signals);
            
            % 检测4个信道中最有可能存在preamble的信道（FFT幅值最大的）
            obj = obj.detectActiveChannel();
            
            % 获得确定信道的信号
            obj.preambleSignal = obj.splitSignal(obj.preambleActiveChannel, :);
            
            % 检测preamble，确定存在preamble并且获得第一个preamble出现的窗口和preamble数目
            obj = obj.detectPreambleBinBehind();
            
            % 通过preamble和SFD的bin来计算CFO和winoffset
            obj = obj.getcfoWinoff();
            
            % 调整信号的winoffset
            signals = circshift(signals, -round(obj.winOffset));
            obj.preambleSignal = circshift(obj.preambleSignal, -round(obj.winOffset));
            
            
%             figure(3);
%             stft(signals(1:60*obj.loraSet.dine), obj.loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',obj.loraSet.fft_x);
            
            % 根据cfo重新生成带有decfo的idealchirp，用于解调
            obj = obj.rebuildIdealchirpCfo(0);
            
%             FFT_plot(obj.preambleSignal(1:16*obj.loraSet.dine), obj.loraSet, obj.cfoDownchirp, 16);
%             FFT_plot(obj.preambleSignal(1:16*obj.loraSet.dine), obj.loraSet, obj.cfoUpchirp, 16);

            obj = obj.getSFDPos();
            
            % 获得payload阶段的跳信道矩阵(TODO: 不一定在窗口的14.25的位置处)
            obj = obj.getDownchirpSync();
            
            % 将信号四个信道滤波后保存在四个流里面
            obj = obj.divideChannel(signals);
            
            % 根据第一个跳信道subchirp（bin为0）来重新对齐
            obj = obj.alignWindows();
            
            % demodulate 解调信号
            obj = obj.demodulateSubchirp();
%             obj = obj.demodulateTest();
        end

        function obj = decodeVote(obj, signals)
            
            % 将信号拆分成channelNum个信道信号
            obj = obj.divideChannel(signals);
            
            % 检测4个信道中最有可能存在preamble的信道（FFT幅值最大的）
            obj = obj.detectActiveChannel();
            
            % 获得确定信道的信号
            obj.preambleSignal = obj.splitSignal(obj.preambleActiveChannel, :);
            
            % 检测preamble，确定存在preamble并且获得第一个preamble出现的窗口和preamble数目
            obj = obj.detectPreambleBinBehind();
            
            % 通过preamble和SFD的bin来计算CFO和winoffset
            obj = obj.getcfoWinoff();
            
            % 调整信号的winoffset
            signals = circshift(signals, -round(obj.winOffset));
            obj.preambleSignal = circshift(obj.preambleSignal, -round(obj.winOffset));
                  
            % 根据cfo重新生成带有decfo的idealchirp，用于解调
            obj = obj.rebuildIdealchirpCfo(0);

            obj = obj.getSFDPos();
            
            % 获得payload阶段的跳信道矩阵(TODO: 不一定在窗口的14.25的位置处)
            obj = obj.getDownchirpSync();
            
            % 将信号四个信道滤波后保存在四个流里面
            obj = obj.divideChannel(signals);
            
            % 根据第一个跳信道subchirp（bin为0）来重新对齐
            obj = obj.alignWindows();
            
            % demodulate 解调信号
%             obj = obj.demodulateSubchirp();
            obj = obj.demodulateVote();
        end

        function output_matrix = generateArray(obj, num, range)
            % 参数
            mask = uint64(2^31);  
            multiplier = uint64(1103515245); 
            c = uint64(12345);  
            
            % 确保输入在0到1023之间
            x = uint64(mod(num, obj.loraSet.fft_x));
            
            % 生成1x1000的矩阵，元素在0到15之间
            output_matrix = uint64(zeros(1, 1000));
            for i = 1:1000
                x = mod(multiplier * x + c, mask);
                tmp = double(x);
                tmp = tmp / 10000;
                tmp = mod(fix(tmp), 10);
        
                output_matrix(i) = uint64(mod(tmp, range) + 1);
            end
        end
        
        function obj = divideChannel(obj, signals)
            obj.splitSignal = [];
            for channel = 1 : obj.channelNum
                signalOut = obj.signalFrequencyShift(signals, -obj.channelList(channel));
                signalOut = obj.lowPassFilterFir(signalOut);
                obj.splitSignal = [obj.splitSignal; signalOut];
            end
        end
        
        function obj = detectActiveChannel(obj)
            dine = obj.loraSet.dine;
            fft_x = obj.loraSet.fft_x;
            max_amp = 0;
            max_bin = 0;
            active_channel = 0;
            for channel = 1:obj.channelNum
                signal_tmp = obj.splitSignal(channel, 1*dine+1:2*dine);
                dechirp = signal_tmp .* obj.idealDownchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine);
                [amp, bin] = max(dechirp_fft);
                if amp > max_amp
                    max_amp = amp;
                    max_bin = bin;
                    active_channel = channel;
                end
            end
            obj.preambleActiveChannel = active_channel;
        end
        
        function obj = detectPreambleBin(obj)
            dine = obj.loraSet.dine;
            fft_x = obj.loraSet.fft_x;
            preamble_len = obj.loraSet.Preamble_length;
            candidate = zeros(1, preamble_len + 2);
            for t = 1:preamble_len + 2
                signal_tmp = obj.preambleSignal((t-1)*dine+1 : t*dine);
                dechirp = signal_tmp .* obj.idealDownchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine);
                [~, candidate(t)] = max(dechirp_fft);
            end
            Preamble_bin = mode(candidate);
            Preamble_start_pos = find(candidate == Preamble_bin);
            Preamble_start_pos = Preamble_start_pos(1);
            Preamble_num = 0;
            for t = Preamble_start_pos:preamble_len + 2
                if candidate(t) == Preamble_bin
                    Preamble_num = Preamble_num + 1;
                else
                    break;
                end
            end
            obj.preambleStartPos = Preamble_start_pos;
            obj.preambleNum = Preamble_num;
        end
        
        % 通过检测 preamble 的众数来确定 preamble 的位置
        function obj = detectPreambleBinBehind(obj)
            obj.DebugUtil.info("\t", "找到最后一个 preamble 的窗口位置");
            dine = obj.loraSet.dine;
            fft_x = obj.loraSet.fft_x;
            preamble_len = obj.loraSet.Preamble_length;
            candidate = zeros(1, preamble_len + 2);

            for t = 1:preamble_len + 4
                signal_tmp = obj.preambleSignal((t-1)*dine+1 : t*dine);
                dechirp = signal_tmp .* obj.idealDownchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine);
                [~, candidate(t)] = max(dechirp_fft);
            end
            % e.g. candidate = [308, 823, 823, 823, 823, 823, 823, 823, 823, 831, 839, 839]
            Preamble_bin = mode(candidate);           % 找到存在cfo 的 preamble 的 bin
            % 找到 sync word 前一个preamble
            for t = 2 : preamble_len + 4
                % 找到符合 sync word 特性的最后一个preamble, syncword 与 preamble 之间的间隔为 8 个 bin
                if (candidate(t) - candidate(t-1) >= 7 && candidate(t) - candidate(t-1) <= 9) ...
                        && abs(Preamble_bin - candidate(t-1)) <= 1
                    Preamble_start_pos = t - 1;
                    break;
                end
            end
            % 处理找不到最后一个preamble的情况
            if exist('Preamble_start_pos', 'var') == 0
                Preamble_start_pos = find(candidate == Preamble_bin);
                Preamble_start_pos = Preamble_start_pos(end);
            end
            obj.preambleEndPos = Preamble_start_pos;   % e.g. 9
            obj.PreambleBin = Preamble_bin;   % e.g. 823
        end
        
        % 利用 Preamble(Base-upchirp) 和 SFD(Base-downchirp) 相反偏移的性质，计算 CFO 和 winoffset
        function obj = getcfoWinoff(obj)
            % 计算主峰的 CFO (需要补零操作, 为了更好地评估峰值频率，可以通过用零填充原始信号来增加分析窗的长度。这种方法以更精确的频率分辨率自动对信号的傅里叶变换进行插值)
            % 对 Preamble 阶段的 FFT 峰值进行排序，得到前 filter 的峰
            zeropadding_size = obj.loraSet.factor;                   % 设置补零的数量，这里的 decim 表示，补上 decim-1 倍窗口的零，计算 FFT 时一共是 decim 倍的窗口（decim+1）, e.g. 16
            d_sf = obj.loraSet.sf;
            d_bw = obj.loraSet.bw;
            dine = obj.loraSet.dine;
            fft_x = obj.loraSet.fft_x;
            Preamble_num = obj.loraSet.Preamble_length;
            downchirp = obj.idealDownchirp;
            upchirp = obj.idealUpchirp;
            filter_num = obj.loraSet.filter_num;
            leakage_width1 = obj.loraSet.leakage_width1;    % 0.0050
            leakage_width2 = obj.loraSet.leakage_width2;    % 0.9950
            signal = obj.preambleSignal;
            preambleEndPosTemp = obj.preambleEndPos;

            dine_zeropadding = dine * zeropadding_size * 2 ^ (10 - d_sf);   % e.g. 16384 * 16 * 2 ^ (10 - 10) = 262144
            fft_x_zeropadding = fft_x * zeropadding_size * 2 ^ (10 - d_sf);  % e.g. 1024 * 16 * 2 ^ (10 - 10) = 16384

            % 获取最后一个preamble窗口的若干个峰值，找到最接近preambleBin的峰
            samples = reshape(signal((preambleEndPosTemp - 8) * dine + 1 : preambleEndPosTemp * dine), [dine, Preamble_num]).';  % e.g. 8 * 16384[]
            samples_fft = abs(fft(samples .* downchirp, dine_zeropadding, 2));  %  e.g. 8 * 262144[]
            samples_fft_merge = samples_fft(:, 1 : fft_x_zeropadding) + samples_fft(:, dine_zeropadding - fft_x_zeropadding + 1 : dine_zeropadding);  % e.g. 8 * 16384[]
            [peak, pos] = sort(samples_fft_merge(1 : Preamble_num, :), 2, 'descend');         % 对 FFT 进行排序
            Peak_pos = zeros(size(pos, 1), filter_num);
            Peak_amp = zeros(size(peak, 1), filter_num);
            Peak_pos(:, 1) = pos(:, 1);
            Peak_amp(:, 1) = peak(:, 1);
            for row = 1 : size(pos, 1)
                temp_array = ones(1, size(pos, 2));
                for list = 1 : filter_num
                    temp_array = temp_array & (abs(Peak_pos(row, list) - pos(row, :)) > fft_x_zeropadding * leakage_width1 & abs(Peak_pos(row, list) - pos(row, :)) < fft_x_zeropadding * leakage_width2);
                    temp_num = find(temp_array == 1, 1, 'first');
                    Peak_pos(row, list + 1) = pos(row, temp_num);
                    Peak_amp(row, list + 1) = peak(row, temp_num);
                end
            end

            % 寻找与第一个窗口的峰（默认第一个窗口只有包1的峰）相近的峰，得到与其相近且重复次数最多的 bin，记作 Preamble 的 bin
            if Peak_pos(2) == Peak_pos(1)
                upchirp_ref = Peak_pos(1);
            else
                upchirp_ref = Peak_pos(2);
            end
            upchirp_index = abs(Peak_pos-upchirp_ref) < fft_x_zeropadding*leakage_width1 | abs(Peak_pos-upchirp_ref) > fft_x_zeropadding*leakage_width2;
            upchirp_bin = (Peak_pos(upchirp_index));
            upchirp_peak = mode(upchirp_bin);

            % 已知 SFD downchirp 的位置，得到 SFD downchirp 的 bin
            SFD_samples = signal((preambleEndPosTemp + 3) * dine + 1 : (preambleEndPosTemp + 4) * dine);
            SFD_samples_fft = abs(fft(SFD_samples .* upchirp, dine_zeropadding));
            samples_fft_merge = SFD_samples_fft(1 : fft_x_zeropadding) + SFD_samples_fft(dine_zeropadding - fft_x_zeropadding + 1 : dine_zeropadding);
            [~, downchirp_peak] = max(samples_fft_merge);   % e.g. 2352
%             figure(1);
%             FFT_plot(signal((preambleEndPosTemp-8)*dine+1:(preambleEndPosTemp+4)*dine), obj.loraSet, downchirp, 12);
%             figure(2);
%             FFT_plot(signal((preambleEndPosTemp-8)*dine+1:(preambleEndPosTemp+4)*dine), obj.loraSet, upchirp, 12);

            % 计算 CFO 和窗口偏移量 (CFO = 2^SF - (preamble bin值 + SFD bin值), 无载波频率偏移时，preamble bin值 + SFD bin值 等于 2^SF)
            if upchirp_peak + downchirp_peak < fft_x_zeropadding*0.5
                cfo_bin = upchirp_peak + downchirp_peak - 2;
                obj.cfo = -cfo_bin/2/fft_x_zeropadding * d_bw;
                obj.winOffset = (downchirp_peak - upchirp_peak) / 2^(11-d_sf);
            elseif upchirp_peak + downchirp_peak > fft_x_zeropadding*1.5
                cfo_bin = upchirp_peak + downchirp_peak - fft_x_zeropadding*2 - 2;
                obj.cfo = -cfo_bin/2/fft_x_zeropadding * d_bw;
                obj.winOffset = (downchirp_peak - upchirp_peak) / 2^(11-d_sf);
            else  % e.g. 13150 + 2352 = 15502
                cfo_bin = upchirp_peak + downchirp_peak - fft_x_zeropadding - 2; % e.g. 15502 - 16384 - 2 = -884
                obj.cfo = -cfo_bin / 2 / fft_x_zeropadding * d_bw;  % e.g. -884 / 2 / 16384 * 125000 = -3.3722e3
                obj.winOffset = (fft_x_zeropadding - (upchirp_peak - downchirp_peak)) / 2 ^ (11 - d_sf);  % e.g. (16384 - (13150 - 2352)) / 2^(11-10) = 5586 / 2 = 2793
            end
            obj.upchirpbin = upchirp_peak;     % e.g. [13150, 13150, 13150, 13150, 13150, 13149, 13149, 13149]
            obj.downchirpbin = downchirp_peak;
        end
        
        function obj = getSFDPos(obj)
            fftX = obj.loraSet.fft_x;
            dine = obj.loraSet.dine;
            preamble_len = obj.loraSet.Preamble_length;
            binRecord = zeros(1, preamble_len + 4);

%             figure(1);
%             FFT_plot(obj.preambleSignal(1:(preamble_len+4)*dine), obj.loraSet, obj.cfoDownchirp, 12);

            for windows = 1 : preamble_len + 4
                signal = obj.preambleSignal((windows-1)*dine+1:windows*dine);
                dechirp = signal .* obj.cfoDownchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
                [~, binRecord(windows)] = max(dechirp_fft);
            end
            for windows = 3 : preamble_len + 4
                first = binRecord(windows-2);
                second = binRecord(windows-1);
                third = binRecord(windows);
                if  (second - first >= 7 && second - first <= 9) && (third - second >= 7 && third - second <= 9)
                    SFDPos = windows;
                    break;
                end
            end
            % 处理找不到最后一个preamble的情况
            if exist('SFDPos', 'var') == 0
               SFDPos = 10;
            end
            obj.SFDPos = SFDPos;
        end
        
        function obj = getDownchirpSync(obj)
            fftX = obj.loraSet.fft_x;
            dine = obj.loraSet.dine;
                        
            signal = obj.preambleSignal((obj.SFDPos+4.25)*dine+1:(obj.SFDPos+5.25)*dine);
            dechirp = signal .* obj.cfoUpchirp;
            dechirp_fft = abs(fft(dechirp, dine));
            dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
            [~, downchirp_bin] = max(dechirp_fft);

%             load("E:\CHchirp\Result\random_record.mat");
            obj.channelArray = obj.channelMatrix(downchirp_bin, :);
        end
        
        function obj = alignWindows(obj)
            fftX = obj.loraSet.fft_x * obj.loraSet.factor;
            dine = obj.loraSet.dine;
            subchirpNum = obj.subchirpNum;
            dechirp_bin_tmp = zeros(1, subchirpNum);
            signal = obj.splitSignal(:, (obj.SFDPos+5.25)*dine+1:(obj.SFDPos+6.25)*dine);
            
            for subchirpCount = 0:subchirpNum-1
                chirpIntegrated = signal(obj.channelArray(subchirpCount+1), subchirpCount*dine/subchirpNum+1 : (subchirpCount+1)*dine/subchirpNum);
                dechirp = obj.cfoDownchirp(subchirpCount*dine/subchirpNum+1 : (subchirpCount+1)*dine/subchirpNum) .* chirpIntegrated;
                dechirp_fft = abs(fft(dechirp, dine*obj.loraSet.factor));
                dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(end-fftX+1:end);
                [~, dechirp_bin_tmp(subchirpCount+1)] = max(dechirp_fft);
            end
            % obj.timeOffset = round(mean(dechirp_bin_tmp));
            obj.timeOffset = dechirp_bin_tmp(end);
%             obj.timeOffset = 656;
        end
        
        function obj = demodulateSubchirp(obj)
            chirpNum = obj.loraSet.payloadNum;
            subchirpNum = obj.subchirpNum;
            fftX = obj.loraSet.fft_x;
            dine = obj.loraSet.dine;
            signal = obj.splitSignal(:, (obj.SFDPos+5.25)*dine-obj.timeOffset+1:end);

%             stft(signal(1, 1:10*obj.loraSet.dine), obj.loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',obj.loraSet.fft_x);

            subchirp_bin = zeros(chirpNum, subchirpNum);
            for chirp_count = 0:chirpNum-1
                dechirp_bin_tmp = zeros(1, subchirpNum);
                for subchirpCount = 0:subchirpNum-1
                    chirpIntegrated = signal(obj.channelArray(chirp_count*subchirpNum+subchirpCount+1), ...
                        chirp_count*dine + subchirpCount*dine/subchirpNum+1 : chirp_count*dine + (subchirpCount+1)*dine/subchirpNum);
                    dechirp = obj.cfoDownchirp(subchirpCount*dine/subchirpNum+1 : (subchirpCount+1)*dine/subchirpNum) .* chirpIntegrated;
                    dechirp_fft = abs(fft(dechirp, dine));
                    dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
                    [~, dechirp_bin_tmp(subchirpCount+1)] = max(dechirp_fft);
                end
                subchirp_bin(chirp_count+1, :) = dechirp_bin_tmp;
            end
            obj.payloadBin = subchirp_bin;
        end

        %{
            投票机制的解subchirp bin
        %}
        function obj = demodulateVote(obj)
            chirpNum = obj.loraSet.payloadNum;
            subchirpNum = obj.subchirpNum;
            factor = obj.loraSet.factor;
            fftX = obj.loraSet.fft_x;
            dine = obj.loraSet.dine;
            % 定位到第一个chirp的窗口
            offset = 2*factor;  % 由于需要搜索窗口，需要给信号前面补充一些信号元素
            splitSignal = obj.splitSignal(:, (obj.SFDPos+5.25)*dine-obj.timeOffset+1 - offset:end);

            % 从第一个窗口开始对后面每一个窗口的subchirp进行解调
            subchirp_bin = zeros(chirpNum, 1);
            for chirp_count = 0:chirpNum-1
                binTmp = zeros(subchirpNum, factor+1);
                % 对每个chirp窗口进行滑动获取bin值
                for slideWindow = -factor/2:factor/2
%                     binTmp = zeros(1, subchirpNum);
                    for subchirpCount = 0:subchirpNum-1
                        chirpIntegrated = splitSignal(obj.channelArray(chirp_count*subchirpNum+subchirpCount+1), ...
                            chirp_count*dine + subchirpCount*dine/subchirpNum+1 + offset + slideWindow : ... 
                            chirp_count*dine + (subchirpCount+1)*dine/subchirpNum + offset + slideWindow);
                        dechirp = obj.cfoDownchirp(subchirpCount*dine/subchirpNum+1 : (subchirpCount+1)*dine/subchirpNum) .* chirpIntegrated;
                        dechirp_fft = abs(fft(dechirp, dine));
                        dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
                        % 获取这个窗口内最大的bin
                        [~, binTmp(subchirpCount+1, slideWindow+factor/2+1)] = max(dechirp_fft);
                    end
                end
                % 投票选出峰值
                subchirp_bin(chirp_count + 1) = mode(binTmp(:));
            end
            obj.payloadBin = subchirp_bin;
        end

    end
end