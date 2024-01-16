classdef NogChirpDecoder < LoraDecoder
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
        upchirpbin;
        downchirpbin;
        fftBinDeWinRecord;   % 记录整个解码窗口的峰值位置
    end

    methods
        % 初始化方法
        function obj = NogChirpDecoder(loraSet)
            obj@LoraDecoder(loraSet);
        end

        function obj = clear(obj)
            % 清空类中的某些中间值
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
            obj.fftBinDeWinRecord = [];
        end

        function obj = decode(obj, signals)
            obj = obj.clear();

            obj.preambleSignal = signals;

            % 检测 preamble，确定存在 preamble 并且获得最后一个 preamble 出现的窗口和 preamble 数目
            obj = obj.detectPreambleBinBehind();

            % 通过 preamble 和 SFD 的 bin 来计算 CFO 和 winoffset
            obj = obj.getcfoWinoff();

            % 调整信号的 winoffset
            obj.preambleSignal = circshift(obj.preambleSignal, -round(obj.winOffset));

            % 根据 cfo 重新生成带有 decfo 的 idealchirp，用于解调
            obj = obj.rebuildIdealchirpCfo(0);

            obj = obj.getSFDPos();

            obj = obj.NogSingleDecode();
            % obj = obj.singleDecode();
        end

        % 方法: payload de-chirp (Default)
        % 参数:
        % 结果: obj.payloadBin
        function obj = singleDecode(obj)
            fftX = obj.loraSet.fft_x;
            dine = obj.loraSet.dine;
            preambleSignalTemp = obj.preambleSignal;
            binArr = zeros(1, obj.loraSet.payloadNum);

            % dechirp做FFT，获得最大峰值即为信道矩阵的信息
            signal = preambleSignalTemp((obj.SFDPos + 2.25)*dine + 1 : end);
            for window_i = 1 : obj.loraSet.payloadNum
                signals = signal((window_i - 1) * dine + 1 : window_i * dine);
                dechirp = signals .* obj.cfoDownchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1 : fftX) + dechirp_fft(dine - fftX + 1 : dine);
                [~, binArr(window_i)] = max(dechirp_fft);
            end

            obj.payloadBin = binArr;
        end

        % 方法: chirp de-chirp (Default)
        % 参数:
        % 结果: [signalOut]
        function [signalOut] = decodeChirp(obj, chirp)
            S_t = obj.cfoDownchirp .* chirp;
            dechirp_fft = abs(fft(S_t, obj.loraSet.dine)); % fft()函数的结果是个复数,取绝对值表示幅值
            signalOut = dechirp_fft(1 : obj.loraSet.fft_x) + dechirp_fft(obj.loraSet.dine - obj.loraSet.fft_x + 1 : obj.loraSet.dine);
        end

        % 方法: 高通滤波, 保留 -bw/2 到 bw/2 的信号
        % 参数:
        % -- chirp: 信号
        % 结果: [signalOutRevise]
        function [signalOutRevise] = filterOutOtherCH(obj, chirp)
            f_pass = obj.loraSet.bw / 2;    % 低通滤波器频率范围, 只保留 -bw/2 到 bw/2 的信号。
            w_pass = f_pass / (obj.loraSet.sample_rate/2);  % 计算归一化频率
            order = 200;   % 滤波器阶数
            b = fir1(order, w_pass, 'low');   % FIR 低通滤波器
            delay = mean(grpdelay(b, obj.loraSet.dine, obj.loraSet.sample_rate));   % 计算滤波器的延迟
            signalIn = chirp;
            signalIn(obj.loraSet.dine + 1 : obj.loraSet.dine + delay) = 0;  % 延迟补零 e.g. dine:(16384 -> 16884)
            signalOutRevise = filter(b, 1, signalIn);  % 滤波

            signalOutRevise(1 : delay) = [];  % 去除前面部分的延迟偏移
        end

        % 方法: 提取候选峰, 筛选能量大于阈值的峰值对应的 bin 值
        % 参数:
        % -- signalInfo: 候选峰信息(峰值, Bin 值)
        % -- winSize: 窗口大小(0 ~ 1)
        % 结果: [signalOutRevise]
        function [groupInfo] = powerExtraction(obj, signalInfo, winSize)
            peak = signalInfo(1, :);
            pos = signalInfo(2, :);
            PowerAmp = obj.loraSet.fft_x;
            PowerThreshold = PowerAmp * winSize / 2; % 阈值为窗口大小能量的一半
            PeakLen = length(signalInfo(1, :));
            for i = 1 : PeakLen
                if peak(1, i) < PowerThreshold
                    peak = peak(1 : i - 1);
                    pos = pos(1 : i - 1);
                    break;
                end
            end
            groupInfo = [peak; pos];
        end

        % 方法: 生成下行扩频信号滑动窗口, 自左向右, winSize 为窗口大小(0 ~ 1), offset 为窗口偏移(0 ~ 2^SF - 1)
        % 参数:
        % -- winSize: 窗口大小(0 ~ 1)
        % -- step: 窗口偏移(0 ~ 2^SF - 1)}
        % 结果: [downchirpSW]
        function [downchirpSW] = downchirpSW(obj, winSize, offset)
            function [noneArr] = noneArr(seg)
                if length(seg) == 1
                    noneArr = [];
                else
                    noneArr = zeros(size(seg));
                end
            end

            cmx = 1 + 1 * 1i;
            pre_dir = 2 * pi;
            d_symbols_per_second = obj.loraSet.bw / obj.loraSet.fft_x;
            T = -0.5 * obj.loraSet.bw * d_symbols_per_second;
            d_samples_per_second = obj.loraSet.sample_rate;       % sdr-rtl的采样率
            d_dt = 1/d_samples_per_second;         % 采样点间间隔的时间
            f0 = obj.loraSet.bw / 2 + obj.cfo;

            r_t = offset / obj.loraSet.fft_x;  % e.g. 1/8 (128/1024)
            seg_1 = d_dt * (0 : 1 : floor(r_t * obj.loraSet.dine)- 1);   % floor(x): x 向下取整
            seg_2 = d_dt * (floor(r_t * obj.loraSet.dine) : 1 : floor((r_t + winSize) * obj.loraSet.dine)- 1);
            seg_3 = d_dt * (floor((r_t + winSize) * obj.loraSet.dine) : 1 : obj.loraSet.dine- 1);
            % seg_1 = d_dt * (0 : 1 : r_t * (obj.loraSet.dine- 1));
            % seg_2 = d_dt * (r_t * (obj.loraSet.dine- 1) : 1 : (r_t + winSize) * (obj.loraSet.dine- 1));
            % seg_3 = d_dt * ((r_t + winSize) * (obj.loraSet.dine- 1) : 1 : obj.loraSet.dine- 1);
            Downchirp_1 =  noneArr(seg_1);
            Downchirp_2 = cmx * (cos(pre_dir .* seg_2 .* (f0 + T * seg_2)) + sin(pre_dir .* seg_2 .* (f0 + T * seg_2)) * 1i);
            Downchirp_3 =  noneArr(seg_3);
            downchirpSW = [Downchirp_1, Downchirp_2, Downchirp_3];
        end

        % 方法: 滑动窗口的方法解调, 自左向右移动 (downchirp 切分)
        % 参数:
        % -- chirp: 信号
        % -- winSize: 窗口大小(0 ~ 1)
        % -- step: 窗口偏移(0 ~ 2^SF - 1)}
        % 结果: obj.fftBinDeWinRecord
        function obj = decodeSildWindow(obj, chirp, winSize, step)
            dechirp_fft = obj.decodeChirp(chirp);
            signalWhole = obj.findpeaksWithShift(dechirp_fft, obj.loraSet.fft_x); % 找峰值
            deWinpeakposInfo = obj.powerExtraction(signalWhole, 1);
            obj.fftBinDeWinRecord = deWinpeakposInfo(2,:);

            % zeropadding_size = obj.loraSet.factor;
            % dine_zeropadding = obj.loraSet.dine * zeropadding_size * 2 ^ (10 - obj.loraSet.sf);   % e.g. 16384 * 16 * 2 ^ (10 - 10) = 262144
            % fft_x_zeropadding = obj.loraSet.fft_x * zeropadding_size * 2 ^ (10 - obj.loraSet.sf);  % e.g. 1024 * 16 * 2 ^ (10 - 10) = 16384

            for i = 0 : 1: ((obj.loraSet.fft_x - (obj.loraSet.fft_x * winSize)) / step)
                % samples = chirp;
                % samples_fft = abs(fft(samples .* obj.downchirpSW(winSize, i * step), dine_zeropadding));  %  e.g. 8 * 262144[]
                % samples_fft_merge = samples_fft(:, 1 : fft_x_zeropadding) + samples_fft(:, dine_zeropadding - fft_x_zeropadding + 1 : dine_zeropadding);  % e.g. 8 * 16384[]
                % samplesOut = obj.findpeaksWithShift(samples_fft_merge, obj.loraSet.fft_x);     % 找峰值
                % [peak, pos] = sort(samples_fft_merge, 'descend');         % 对 FFT 进行排序

                S_t = obj.downchirpSW(winSize, i * step) .* chirp;
                dechirp_fft = abs(fft(S_t, obj.loraSet.dine));      % fft()函数的结果是个复数, 取绝对值表示幅值
                dechirp_fft = dechirp_fft(1 : obj.loraSet.fft_x) + dechirp_fft(obj.loraSet.dine - obj.loraSet.fft_x + 1 : obj.loraSet.dine);
                sildWinOut = obj.findpeaksWithShift(dechirp_fft, obj.loraSet.fft_x); % 找峰值
                sildWinPEOut = obj.powerExtraction(sildWinOut, winSize);
                % [sildWinPeak] = sildWinPEOut(1, :);
                [sildWinPos] = sildWinPEOut(2, :);

                % disp(['sildWinPos-', num2str(i), '-: [', num2str(sildWinPos), ']']);

                recordLen = length(obj.fftBinDeWinRecord);
                % 提取候选峰, 筛选能量大于阈值的峰值对应的 bin 值
                for j = 1 : recordLen
                    element = obj.fftBinDeWinRecord(j);
                    if ~ismember(element, sildWinPos) && ~ismember(element + 1, sildWinPos) && ~ismember(element - 1, sildWinPos)  % e.g. fftBinDeWinRecord = [810(目标bin), 555(噪声), 467(非正交bin)]  pos = [810(目标bin), 244(噪声)]
                        index_to_remove = (obj.fftBinDeWinRecord == element);
                        obj.fftBinDeWinRecord(index_to_remove) = -1;
                    end
                end
            end
            % 找到非负值的元素的索引
            non_negative_indices = (obj.fftBinDeWinRecord >= 0);
            obj.fftBinDeWinRecord = obj.fftBinDeWinRecord(non_negative_indices);
        end

        % 方法: 滑动窗口的方法解调, 自左向右移动 (datachirp 切分)
        % 参数:
        % -- chirp: 信号
        % -- winSize: 窗口大小(0 ~ 1)
        % -- step: 窗口偏移(0 ~ 2^SF - 1)}
        % 结果: obj.fftBinDeWinRecord
        function obj = decodeChirpSildWin(obj, chirp, winSize, step)
            dechirp_fft = obj.decodeChirp(chirp);
            signalWhole = obj.findpeaksWithShift(dechirp_fft, obj.loraSet.fft_x); % 找峰值
            deWinpeakposInfo = obj.powerExtraction(signalWhole, 1);
            obj.fftBinDeWinRecord = deWinpeakposInfo(2,:);

            for i = 0 : 1: ((obj.loraSet.fft_x - (obj.loraSet.fft_x * winSize)) / step)
                chirpTmp = chirp;
                chirpStart = (i * step / obj.loraSet.fft_x) * obj.loraSet.dine + 1;  % 确定 chirp 切片的起始位置
                chirpEnd = chirpStart + winSize * obj.loraSet.dine;  % 确定 chirp 切片的结束位置
                chirpTmp(1 : chirpStart - 1) = 0;               % 将 chirp 切片之前的位置置零
                chirpTmp(chirpEnd : obj.loraSet.dine) = 0;      % 将 chirp 切片之后的位置置零

                S_t = obj.cfoDownchirp .* chirpTmp;
                dechirp_fft = abs(fft(S_t, obj.loraSet.dine));      % fft()函数的结果是个复数, 取绝对值表示幅值
                dechirp_fft = dechirp_fft(1 : obj.loraSet.fft_x) + dechirp_fft(obj.loraSet.dine - obj.loraSet.fft_x + 1 : obj.loraSet.dine);
                sildWinOut = obj.findpeaksWithShift(dechirp_fft, obj.loraSet.fft_x); % 找峰值
                sildWinPEOut = obj.powerExtraction(sildWinOut, winSize);
                % [sildWinPeak] = sildWinPEOut(1, :);
                [sildWinPos] = sildWinPEOut(2, :);

                % disp(['sildWinPos-', num2str(i), '-: [', num2str(sildWinPos), ']']);

                recordLen = length(obj.fftBinDeWinRecord);
                % 提取候选峰, 筛选能量大于阈值的峰值对应的 bin 值
                for j = 1 : recordLen
                    element = obj.fftBinDeWinRecord(j);
                    if ~ismember(element, sildWinPos)  % e.g. fftBinDeWinRecord = [810(目标bin), 555(噪声), 467(非正交bin)]  pos = [810(目标bin), 244(噪声)]
                        index_to_remove = (obj.fftBinDeWinRecord == element);
                        obj.fftBinDeWinRecord(index_to_remove) = -1;
                    end
                end
            end
            % 找到非负值的元素的索引
            non_negative_indices = (obj.fftBinDeWinRecord >= 0);
            obj.fftBinDeWinRecord = obj.fftBinDeWinRecord(non_negative_indices);
        end

        % 方法: 非正交信道解码
        % 参数:
        % 结果: obj.payloadBin
        function obj = NogSingleDecode(obj)
            dine = obj.loraSet.dine;
            preambleSignalTemp = obj.preambleSignal;
            % binArr = zeros(1, obj.loraSet.payloadNum);

            % dechirp做FFT，获得最大峰值即为信道矩阵的信息
            signal = preambleSignalTemp((obj.SFDPos + 2.25)*dine + 1 : end);
            for window_i = 5 : obj.loraSet.payloadNum
                windowChirp = signal((window_i - 1) * dine + 1 : window_i * dine);
                
                save('TwoCollision_Chirp_1.mat',"windowChirp");
                break;
                signalOutRevise = obj.filterOutOtherCH(windowChirp);
                obj = obj.decodeSildWindow(signalOutRevise, 1/4, 128);   % 信号, 窗口大小, 步长
                % obj = obj.decodeChirpSildWin(signalOutRevise, 1/8, 128);   % 信号, 窗口大小, 步长
                % binArr(window_i) = obj.fftBinDeWinRecord;
                disp(obj.fftBinDeWinRecord);
            end
            % obj.payloadBin = binArr;
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

        % 通过 SFD 检测 preamble 的位置
        function obj = detectSFDofCrossCorrelation(obj)

        end

        % 通过检测 preamble 的众数来确定 preamble 的位置
        function obj = detectPreambleBinBehind(obj)
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
            obj.preambleBin = Preamble_bin;   % e.g. 823
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
            binRecordSFD = zeros(1, preamble_len + 4);

%             figure(1);
%             FFT_plot(obj.preambleSignal(1:(preamble_len+4)*dine), obj.loraSet, obj.cfoDownchirp, 12);

            for windows = 1 : preamble_len + 4
                signal = obj.preambleSignal((windows - 1) * dine + 1 : windows * dine);
                dechirp = signal .* obj.cfoDownchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1 : fftX) + dechirp_fft(dine - fftX + 1 : dine);
                [~, binRecordSFD(windows)] = max(dechirp_fft);
            end
            for windows = 3 : preamble_len + 4
                first = binRecordSFD(windows - 2);
                second = binRecordSFD(windows - 1);
                third = binRecordSFD(windows);
                if  (second - first >= 7 && second - first <= 9) && (third - second >= 7 && third - second <= 9)
                    SFDPosTemp = windows;
                    break;
                end
            end
            % 处理找不到最后一个 preamble 的情况
            if exist('SFDPos', 'var') == 0
                SFDPosTemp = 10;
            end
            obj.SFDPos = SFDPosTemp;
        end
    end
end