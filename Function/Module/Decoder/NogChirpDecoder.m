classdef NogChirpDecoder < LoraDecoder
    properties
        payloadBin;
        BinRecord;  % 记录临时的bin值
        preambleSignal;
        preambleStartPos;
        preambleBin;
        preambleEndPos;
        SFDPos;
        window;
        errorFlag;
        errorMsg;
        peakStandard;
        SFDPeakAmp;
        DebugUtil;
        preamblePeak;
        upchirpbin;
        downchirpbin;
        fftBinDeWinRecord;     % 记录整个解码窗口的峰值位置
        detectArray;           % 存放检测到的 bin
        detectArrayCount;      % 计数器
        detectArrayNum;        % 记录队列中存放有效数据的数目
    end

    methods
        % 初始化方法
        function obj = NogChirpDecoder(loraSet)
            obj@LoraDecoder(loraSet);
        end

        function obj = clear(obj)
            % 清空类中的某些中间值
            obj.payloadBin = {};
            obj.BinRecord = [];
            obj.preambleSignal = [];
            obj.preambleStartPos = [];
            obj.preambleEndPos = [];
            obj.SFDPos = [];
            obj.peakStandard = [];
            obj.SFDPeakAmp = [];
            obj.cfo = [];
            obj.winOffset = [];
            obj.cfoDownchirp = [];
            obj.cfoUpchirp = [];
            obj.preamblePeak = [];
            obj.errorFlag = false;  % 每次进行setArgs相当于一次初始化，所以要对errorFlag置为false
            obj.errorMsg = "";
            obj.fftBinDeWinRecord = [];
            obj.detectArray = zeros(1, obj.loraSet.fft_x);   % 存放信道检测到的峰值bin（bin 最多为 fft_x 种）
            obj.detectArrayCount = zeros(1, obj.loraSet.fft_x);  % 计数器
            obj.detectArrayNum = 0;  % 记录信道队列中存放有效数据的数目
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

            obj = obj.OverlapChDecoder();
            % obj = obj.singleDecode();
        end

        %% 方法: 解码测试
        % 参数:
        % 结果:
        function obj = decodeTwoCH(obj, signals)
            obj = obj.clear();

            for Ch_i = 1 : 2
                disp("Channel: -" + num2str(Ch_i) + "-");
                if Ch_i == 2
                    signals = obj.rechangeSignalFreq(signals, - obj.loraSet.bw / 4);  % 中心频率对齐第二信道
                end
                windowsNum = fix(length(signals) / obj.loraSet.dine); % 计算需要移动多少个窗口进行信号检测和解调
                % FFT_plot(obj.preambleSignal, obj.loraSet, obj.idealDownchirp, 23);
                for window_i = 1 : windowsNum % 开始扫描每一个窗口，发现其中的 premable
                    windowChirp = signals((window_i - 1) * obj.loraSet.dine + 1 : window_i * obj.loraSet.dine);
                    [obj, preambleBinTmp] = obj.detect(windowChirp);
                    if ~isempty(preambleBinTmp)
                        disp("♦ preamble bin: [" + regexprep(num2str(preambleBinTmp), '\s+', ' ') + "]");
                        % 根据detect的结果进行解调
                        for num = 1:length(preambleBinTmp)
                            obj.window = window_i;
                            obj.preambleBin = preambleBinTmp(num);
                            obj.preambleSignal = signals;
                            obj = obj.detectPreambleEndPosBehind();  % 检测 preamble 结束位置
                            disp("  • preambleEndPos: " + num2str(obj.preambleEndPos));

                            % begin decoding code
                            obj = obj.getcfoWinoff();   % 通过 preamble 和 SFD 的 bin 来计算 CFO 和 winoffset
                            disp("  • CFO: " + num2str(obj.cfo));
                            disp("  • winOffset: " + num2str(obj.winOffset));
                            if abs(obj.cfo) >= 10e3  % 如果 CFO 太大，代表是另外一个信道的信号，此时直接跳过
                                disp("  • Warning: CFO is too large, skip this signal");
                                fprintf('\n');
                            else
                                obj.preambleSignal = circshift(obj.preambleSignal, -round(obj.winOffset));  % 根据 winoffset 调整信号
                                obj = obj.rebuildIdealchirpCfo(0);  % 根据 cfo 重新生成带有 decfo 的 idealchirp，用于解调
                                obj = obj.getSFDPos();   % 获取 SFD 的位置
                                disp("  • SFD Position: " + num2str(obj.SFDPos));
                                fprintf('\n');
                                obj = obj.OverlapChDecoder(Ch_i);  % 解调
                                obj.payloadBin{end + 1} = obj.BinRecord;
                                obj.BinRecord = [];
                                % end decoding code
                            end
                        end
                    end
                end
            end
        end

        %% 方法: 重叠信道解码
        % 参数:
        % 结果: obj.payloadBin
        function obj = OverlapChDecoder(obj, ChNum)
            dine = obj.loraSet.dine;
            preambleSignalTemp = obj.preambleSignal;
            % binArr = zeros(1, obj.loraSet.payloadNum);

            % dechirp 做 FFT，获得最大峰值即为信道矩阵的信息
            signal = preambleSignalTemp((obj.SFDPos + 2.25) * dine + 1 : end);
            for window_i = 1 : obj.loraSet.payloadNum
                windowChirp = signal((window_i - 1) * dine + 1 : window_i * dine);
                % filename = sprintf('TwoCollision_Chirp_%d.mat', window_i-1);
                % save(filename, 'windowChirp');

                %% 信号正方向 =》 偏移量信号滤波 =》 筛选 Bin 值
                % obj = obj.decodeNonoverlap(windowChirp, ChNum);

                %% 重叠窗口和非重叠窗口解码后 Bin 值交集
                obj = obj.decodeIntersection(windowChirp, ChNum);

                %% 滑动窗口的方法解调
                % signalOutRevise = obj.filterOutOtherCH(windowChirp);
                % obj = obj.decodeSildWindow(signalOutRevise, 1/8, 128);   % 信号, 窗口大小, 步长
                % obj = obj.decodeChirpSildWin(signalOutRevise, 1/8, 128);   % 信号, 窗口大小, 步长
                % binArr(window_i) = obj.fftBinDeWinRecord;
                % disp(obj.fftBinDeWinRecord);
            end
            % obj.payloadBin = binArr;
        end

        %% 方法: payload de-chirp (Default)
        % 参数:
        % 结果: obj.payloadBin
        function obj = singleDecode(obj)
            fft_xTmp = obj.loraSet.fft_x;
            dine = obj.loraSet.dine;
            preambleSignalTemp = obj.preambleSignal;
            binArr = zeros(1, obj.loraSet.payloadNum);

            % dechirp做FFT，获得最大峰值即为信道矩阵的信息
            signal = preambleSignalTemp((obj.SFDPos + 2.25)*dine + 1 : end);
            for window_i = 1 : obj.loraSet.payloadNum
                signals = signal((window_i - 1) * dine + 1 : window_i * dine);
                dechirp = signals .* obj.cfoDownchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1 : fft_xTmp) + dechirp_fft(dine - fft_xTmp + 1 : dine);
                [~, binArr(window_i)] = max(dechirp_fft);
            end

            obj.payloadBin = binArr;
        end

        %% 方法: 生成指定频率偏移的 downchirp
        % 参数:
        % -- Freq: 频率
        % 结果: [DownchirpOut]
        function [DownchirpOut] = cfoFreqShiftDownchirp(obj, Freq)
            cmx = 1 + 1 * 1i;
            pre_dir = 2 * pi;
            Ns = obj.loraSet.bw / obj.loraSet.fft_x;  % 1s 内的 symbol 数 (Ns = 1 / Ts = bw / 2^SF)
            T = -0.5 * obj.loraSet.bw * Ns;
            d_dt = 1 / obj.loraSet.sample_rate;       % 采样点间间隔的时间
            t = d_dt * (0 : 1 : obj.loraSet.dine - 1);
            f0 = obj.loraSet.bw / 2 + obj.cfo + Freq;

            DownchirpOut = cmx * (cos(pre_dir .* t .* (f0 + T * t)) + sin(pre_dir .* t .* (f0 + T * t)) * 1i);
            % r_t = Freq / obj.loraSet.bw;
            % seg_1 = d_dt * (0 : 1 : (1 - r_t) * (obj.loraSet.dine- 1)); % 0~T1, 即 chirp 第 1 段
            % seg_2 = d_dt * ((1 - r_t) * (obj.loraSet.dine- 1) : 1 : obj.loraSet.dine- 1); % T1~1, 即 chirp 第 2 段
            % dataDownchirp_1 = cmx * (cos(pre_dir .* seg_1 .* (f0 + T * seg_1 - Freq)) + sin(pre_dir .* seg_1 .* (f0 + T * seg_1 - Freq)) * 1i);
            % dataDownchirp_2 = cmx * (cos(pre_dir .* seg_2 .* (f0 + T * seg_2 - Freq + obj.loraSet.bw)) + sin(pre_dir .* seg_2 .* (f0 + T * seg_2 - Freq + obj.loraSet.bw)) * 1i);
            % DownchirpOut = [dataDownchirp_1, dataDownchirp_2];
        end

        %% 方法: de-chirp
        % 参数:
        % -- chirp: 信号
        % -- downchirpPick: 下行 chirp 类型选择
        % -- Freq: 频率
        % 结果: [signalOut]
        function [signalOut] = decodeChirp(obj, chirp, downchirpPick, Freq)
            if nargin < 3 || (nargin == 3 && downchirpPick == 1)
                down_chirp = obj.cfoDownchirp;
            elseif nargin == 3 && downchirpPick == 2
                down_chirp = obj.idealDownchirp;
            elseif nargin == 4 && downchirpPick == 3
                down_chirp = obj.cfoFreqShiftDownchirp(Freq);
            end
            dineTmp  = obj.loraSet.dine;
            fft_xTmp = obj.loraSet.fft_x;
            S_t = down_chirp .* chirp;
            dechirp_fft = abs(fft(S_t, dineTmp)); % fft()函数的结果是个复数, 取绝对值表示幅值
            signalOut = dechirp_fft(1 : fft_xTmp) + dechirp_fft(dineTmp - fft_xTmp + 1 : dineTmp);
        end

        %% 方法: 提高/降低 信号整体频率
        % 参数:
        % -- chirp: 信号
        % -- shiftFreq：提高/降低 的频率
        % 结果: [shiftedSignal]
        function [shiftedSignal] = rechangeSignalFreq(obj, chirp, shiftFreq)
            pre_dir = 2 * pi;
            d_dt = 1 / obj.loraSet.sample_rate;
            t = d_dt * (0 : 1: length(chirp) - 1); % 时间序列
            shiftedSignal = chirp .* exp(1i * pre_dir * t * shiftFreq); % 信号整体频率提高/降低 shiftFreq
        end

        %% 方法: 低通/带通/高通 滤波, 保留 satrtFreq 到 endFreq 的信号
        % 参数:
        % -- chirp: 信号
        % -- order: 滤波器阶数
        % -- satrtFreq: 起始频率
        % -- endFreq: 结束频率
        % 结果: [signalOutRevise]
        function [signalOutRevise] = filterOutOtherCH(obj, chirp, order, satrtFreq, endFreq)
            if nargin < 3   % 如果没有提供对应输入参数，则设置默认值
                order = 200;   % 滤波器阶数
                satrtFreq = 0;
                endFreq = obj.loraSet.bw / 2;
            elseif nargin == 3
                satrtFreq = 0;
                endFreq = obj.loraSet.bw / 2;
            end

            f_pass = [satrtFreq endFreq];    % 通带滤波器频率范围, 只保留 satrtFreq 到 endFreq 的信号。
            w_pass = f_pass / (obj.loraSet.sample_rate/2);  % 计算归一化频率

            if satrtFreq ~= 0 && endFreq ~= 0
                b = fir1(order, [w_pass(1) w_pass(2)], 'bandpass');   % FIR 带通滤波器
            elseif satrtFreq == 0 && endFreq ~= 0
                b = fir1(order, w_pass(2), 'low');   % FIR 低通滤波器
            elseif satrtFreq ~= 0 && endFreq == 0
                b = fir1(order, w_pass(1), 'high');   % FIR 高通滤波器
            end

            delay = mean(grpdelay(b, obj.loraSet.dine, obj.loraSet.sample_rate));   % 计算滤波器的延迟
            signalIn = chirp;
            signalIn(obj.loraSet.dine + 1 : obj.loraSet.dine + delay) = 0;  % 延迟补零 e.g. dine:(16384 -> 16884)
            signalOutRevise = filter(b, 1, signalIn);  % 滤波

            signalOutRevise(1 : delay) = [];  % 去除前面部分的延迟偏移
            % signalOutRevise(obj.loraSet.dine - delay + 1 : obj.loraSet.dine) = 1;
            % signalOutRevise = signalOutRevise(1 : obj.loraSet.dine);   % 截取原始信号长度
        end

        %% 方法: 提取候选峰, 筛选能量大于阈值的峰值对应的 bin 值
        % 参数:
        % -- signalInfo: 候选峰信息(峰值, Bin 值)
        % -- winSize: 窗口大小(0 ~ 1)
        % 结果: [groupInfo]
        function [groupInfo] = powerExtraction(obj, signalInfo, winSize)
            peak = signalInfo(1, :);
            pos = signalInfo(2, :);
            PowerAmp = obj.loraSet.fft_x * 2;  % 信号能量幅值
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

        %% 方法: 生成下行扩频信号滑动窗口, 自左向右, winSize 为窗口大小(0 ~ 1), offset 为窗口偏移(0 ~ 2^SF - 1)
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
            Ns = obj.loraSet.bw / obj.loraSet.fft_x;  % 1s 内的 symbol 数 (Ns = 1 / Ts = bw / 2^SF)
            T = -0.5 * obj.loraSet.bw * Ns;
            d_dt = 1 / obj.loraSet.sample_rate;       % 采样点间间隔的时间
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

        %% 方法: 滑动窗口的方法解调, 自左向右移动 (downchirp 切分)
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

        %% 方法: 滑动窗口的方法解调, 自左向右移动 (datachirp 切分)
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

        %% 方法: 信号正方向 =》滤波 =》筛选峰值能量大于原 1/4 的 bin值
        % 参数:
        % -- chirp: 信号
        % -- PowerAmp: 能量幅值
        % -- cfoFlag: 是否带 cfo
        % 结果: obj.payloadBin
        function obj = decodeNonoverlap(obj, chirp, ChNum)
            % if ChNum == 1
            %     FreqShift = -obj.loraSet.bw / 4;
            % elseif ChNum == 2
            %     FreqShift = obj.loraSet.bw / 4;
            % end
            % FreqFliterStart = obj.loraSet.bw / 2;
            % FreqFliterEnd = 0;

            if ChNum == 1
                FreqShift = obj.loraSet.bw / 2;
            elseif ChNum == 2
                FreqShift = -obj.loraSet.bw / 2;
            end
            FreqFliterStart = 0;
            FreqFliterEnd = obj.loraSet.bw / 4;

            dechirp_fft = obj.decodeChirp(chirp);
            signalWhole = obj.findpeaksWithShift(dechirp_fft, obj.loraSet.fft_x);  % 提取候选峰
            deWinpeakposInfo = obj.powerExtraction(signalWhole, 1);
            [OrgsignalPos] = deWinpeakposInfo(2, :);
            [OrgsignalPeak] = deWinpeakposInfo(1, :);
            % disp(['Filter peaks: [', num2str(OrgsignalPeak), ']']);
            % disp(['Candidate pos: [', num2str(OrgsignalPos), ']']);

            chirpFreqIncre = obj.rechangeSignalFreq(chirp, FreqShift);

            chirpFreqIncreRevise = obj.filterOutOtherCH(chirpFreqIncre, 200, FreqFliterStart, FreqFliterEnd);  % 滤波

            dechirp_fft = obj.decodeChirp(chirpFreqIncreRevise, 3, -FreqShift);
            signalInfo = obj.findpeaksWithShift(dechirp_fft, obj.loraSet.fft_x);
            signalPEOut = obj.powerExtraction(signalInfo, 1/4);   % 提取候选峰, 筛选能量大于阈值的峰值对应的 bin 值
            [signalPos] = signalPEOut(2, :);
            [signalPeak] = signalPEOut(1, :);
            % disp(['Filter peaks: [', num2str(signalPeak), ']']);
            % disp(['Filter pos: [', num2str(signalPos), ']']);
            % subplot(211);
            % stft(chirp, obj.loraSet.sample_rate, 'Window', kaiser(64, 2), 'OverlapLength', 32, 'FFTLength', obj.loraSet.fft_x);
            % ylim([-0.1875, 0.1875]);
            % subplot(212);
            % stft(chirpFreqIncreRevise, obj.loraSet.sample_rate, 'Window', kaiser(64, 2), 'OverlapLength', 32, 'FFTLength', obj.loraSet.fft_x);
            % ylim([-0.1875, 0.1875]);
            for i = 1 : length(signalPos)
                pos_i = OrgsignalPos == signalPos(i);
                if any(OrgsignalPeak(pos_i) < signalPeak(i) * 5) && any(OrgsignalPeak(pos_i) > signalPeak(i) * 10 / 3)
                    obj.BinRecord = [obj.BinRecord signalPos(i)];
                end
            end
        end

        %% 重叠窗口和非重叠窗口解码后 Bin 值交集
        % 参数:
        % -- chirp: 信号
        % -- PowerAmp: 能量幅值
        % 结果: obj.BinRecord
        function obj = decodeIntersection(obj, chirp, ChNum)
            if ChNum == 1
                FreqShift = obj.loraSet.bw / 2;
            elseif ChNum == 2
                FreqShift = -obj.loraSet.bw / 2;
            end

            % dechirp_fft = obj.decodeChirp(chirp);
            % signalWhole = obj.findpeaksWithShift(dechirp_fft, obj.loraSet.fft_x);  % 提取候选峰
            % deWinpeakposInfo = obj.powerExtraction(signalWhole, 1);
            % [OrgsignalPos] = deWinpeakposInfo(2, :);
            % [OrgsignalPeak] = deWinpeakposInfo(1, :);
            % disp(['Candidate Bins: [', num2str(OrgsignalPos), ']']);
            % disp(['Corresponding Peaks: [', num2str(OrgsignalPeak), ']']);

            chirpFreqCHG = obj.rechangeSignalFreq(chirp, FreqShift);  % 将信号整体频率变化 FreqShift, 方便后面滤波, eg: (-bw/2, bw/2) --> (0, bw)

            chirpNonOverlapping = obj.filterOutOtherCH(chirpFreqCHG, 150, 0, obj.loraSet.bw / 4);  % 非重叠窗口滤波，eg: (0, bw/4)
            chirpOverlapping = obj.filterOutOtherCH(chirpFreqCHG, 150, obj.loraSet.bw / 4, obj.loraSet.bw); % 重叠窗口滤波，eg: (bw/4, bw)

            dechirp_fft = obj.decodeChirp(chirpNonOverlapping, 3, -FreqShift);  % 非重叠窗口解码, (信号, downchirp 类型, 频率偏移)
            signalInfo = obj.findpeaksWithShift(dechirp_fft, obj.loraSet.fft_x);
            signalPEOut = obj.powerExtraction(signalInfo, 1/4);   % 提取候选峰, 筛选能量大于阈值的峰值对应的 bin 值
            [NonOverlappingPos] = signalPEOut(2, :);
            % [NonOverlappingPeak] = signalPEOut(1, :);

            dechirp_fft = obj.decodeChirp(chirpOverlapping, 3, -FreqShift);  % 重叠窗口解码
            signalInfo = obj.findpeaksWithShift(dechirp_fft, obj.loraSet.fft_x);
            signalPEOut = obj.powerExtraction(signalInfo, 3/4);   % 提取候选峰, 筛选能量大于阈值的峰值对应的 bin 值
            [OverlappingPos] = signalPEOut(2, :);
            % [OverlappingPeak] = signalPEOut(1, :);

            obj.BinRecord = [obj.BinRecord intersect(NonOverlappingPos, OverlappingPos)];   % 重叠窗口和非重叠窗口 Bin 交集
        end

        %% 方法: 检测 preamble 的可能值
        % 参数:
        % -- chirp: 窗口信号
        % 结果: obj.preambleBin
        %% √
        function [obj, preambleBinTmp] = detect(obj, chirp)
            preambleBinTmp = [];
            % preambleBinTmp 用于存放检测出来的lora包的preamble bin
            dine = obj.loraSet.dine;
            fft_x = obj.loraSet.fft_x;

            if isActive(obj, chirp)  % 如果这个信道存在信号
                % 对这个信道的信号乘以downchirp，做FFT，通过findpeaks函数找到突出峰值对应的bin，来判断preamble特性
                chirpTmp = chirp;
                dechirp = chirpTmp .* obj.idealDownchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1 : fft_x) + dechirp_fft(dine - fft_x + 1 : dine);
                % 找出若干个峰值
                result = obj.findpeaksWithShift(dechirp_fft, fft_x);
                deWinpeakposInfo = obj.powerExtraction(result, 1);  % 滤掉低于能量阈值的峰值
                binPos = deWinpeakposInfo(2, :);
                % disp("binPos: [" + regexprep(num2str(binPos), '\s+', ' ') + "]");

                % 对每一个筛选过的峰值进行记录
                detectArrayTmp = obj.detectArray;   % 临时存放信道检测到的峰值bin（bin 最多为 fft_x 种）
                detectAyyayCountTmp = obj.detectArrayCount;   % 临时计数器
                detectArrayNumTmp = obj.detectArrayNum;   % 记录信道队列中存放有效数据的数目

                for pos_i = 1:length(binPos)
                    % disp("pos_i: " + pos_i);
                    arr = obj.detectArray(1 : detectArrayNumTmp);
                    indices = obj.findBin(arr, binPos(pos_i), fft_x);
                    % 如果检测队列为空或者检测队列中没有该元素或没有与该元素绝对值相差1的值，直接加入
                    if detectArrayNumTmp == 0 || isempty(indices)
                        % 检测队列的总数量加1
                        detectArrayNumTmp = detectArrayNumTmp + 1;
                        % 将峰对应的位置加入到检测队列中
                        detectArrayTmp(detectArrayNumTmp) = binPos(pos_i);
                        % 并使其对应的count置1
                        detectAyyayCountTmp(detectArrayNumTmp) = 1;
                    else % 如果检测队列有该元素或与该元素绝对值相差1的值，对其记录的count加1
                        % 找到对应的检测队列中的位置
                        % find_pos = find(obj.detectArray(channel, 1:obj.detectArrayNum(channel)) == binPos(pos_i));
                        % 使其对应的count加1
                        detectAyyayCountTmp(indices) = detectAyyayCountTmp(indices) + 1;
                        % 修改纠正这个检测队列中的bin（调整1个bin）
                        detectArrayTmp(indices) = binPos(pos_i);
                    end
                end
                obj.detectArray = detectArrayTmp;
                obj.detectArrayNum = detectArrayNumTmp;
                obj.detectArrayCount = detectAyyayCountTmp;

                % 将在检测队列中未能出现在此次 preamble 的峰值给剔除掉
                % 对检测队列中的每一个峰值进行遍历，找到没有在此次过滤出的峰出现的
                detectArrayNumTmp = obj.detectArrayNum;
                countTmp = 1;
                for pos_i = 1 : detectArrayNumTmp
                    % 发现检测队里中的峰值在此次中未出现
                    binValue = obj.detectArray(countTmp);
                    indices = obj.findBin(binPos, binValue, fft_x);
                    if isempty(indices)
                        % 在检测队列中进行删除
                        obj.detectArray(countTmp : end) = [obj.detectArray(countTmp + 1 : end), 0];
                        obj.detectArrayCount(countTmp : end) = [obj.detectArrayCount(countTmp + 1 : end), 0];
                        obj.detectArrayNum = obj.detectArrayNum - 1;
                    else
                        countTmp = countTmp + 1;
                    end
                end
                % 将检测队列中 count 数目等于 countThreshold 的元素放入检测成功的队列中
                detectArrayNumTmp = obj.detectArrayNum;
                countTmp = 1;
                countThreshold = 5;
                for pos_i = 1:detectArrayNumTmp
%                     if obj.detectArrayCount(channel, countTmp) == 7
                    if obj.detectArrayCount(countTmp) == countThreshold
                        % 将挑选出的lora包的preamble bin记录
                        preambleBinTmp = [preambleBinTmp, obj.detectArray(countTmp)];
                        % 将detect队列中对应元素删除，和上面做法相同
                        obj.detectArray(countTmp:end) = [obj.detectArray(countTmp+1:end), 0];
                        obj.detectArrayCount(countTmp:end) = [obj.detectArrayCount(countTmp+1:end), 0];
                        obj.detectArrayNum = obj.detectArrayNum - 1;
                    else
                        countTmp = countTmp + 1;
                    end
                end
            end
            % TODO: 如果没有active则清除队列
        end

        %% 方法: preamble detection (未使用)
        % 参数: ~
        % 结果: obj.preambleBin  obj.preambleEndPos
        %% ×
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

        %% 方法: 通过检测 preamble 的众数来确定 preamble 的位置
        % 参数: ~
        % 结果: obj.preambleBin  obj.preambleEndPos
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

        %% 方法: 检测当前信道下对应premableBin的最后一个preamble的位置（为后续cfo计算准备）
        % TODO: 因为涉及到多解码，建议把 preambleEndPos 信息作为输出传出
        % 参数: ~
        % 结果: obj.preambleEndPos
        %% √
        function obj  = detectPreambleEndPosBehind(obj)
            dine = obj.loraSet.dine;
            fft_x = obj.loraSet.fft_x;
            preamble_len = obj.loraSet.Preamble_length;
            windowTmp = obj.window;
            preambleSignalTmp = obj.preambleSignal;
            preambleBinTmp = obj.preambleBin;

            % 用于存储preambleLength+4数目窗口下，每个窗口的峰值
            candidate = cell(1, preamble_len);
            cancidateAmp = cell(1, preamble_len);

            % 搜寻第7个preamble前后的chirp
            for t = windowTmp : windowTmp + 7
                % 每一个窗口做FFT后，将获得的结果找出若干个峰值
                signal_tmp = preambleSignalTmp(t*dine+1 : (t+1)*dine);
                dechirp = signal_tmp .* obj.idealDownchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine);
                % 找出若干个峰值
                result = obj.findpeaksWithShift(dechirp_fft, fft_x);
                % 记录峰值
                candidate{t - windowTmp + 1} = result(2, :);
                cancidateAmp{t - windowTmp + 1} = result(1, :);
            end

            % 找到每个窗口中，最接近preambleBin的值
            [BinArray, AmpArray] = obj.findClosetSyncWordBin(candidate, cancidateAmp, preambleBinTmp);

            % 根据sync word的特征找到最后一个preamble
            Preamble_end_pos = obj.findPosWithSyncWord(BinArray, AmpArray);

            % 处理找不到的情况
            if exist('Preamble_end_pos', 'var') == 0 || Preamble_end_pos == 0
                return;
            else
                obj.preambleEndPos = windowTmp + Preamble_end_pos - 2;
                % 获取preamble的能量
                signal_tmp = preambleSignalTmp((obj.preambleEndPos-1)*dine+1 : (obj.preambleEndPos)*dine);
                dechirp = signal_tmp .* obj.idealDownchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine);
                % 找出若干个峰值，峰值间间隔为leakWidth
                result = obj.findpeaksWithShift(dechirp_fft, fft_x);
                binPos = result(2, :);
                % 记录峰值
                [~, closest_idx] = min(abs(binPos - preambleBinTmp));
                obj.preamblePeak = result(1, closest_idx);
            end
        end

        %% 方法: 从输入的元胞数组中，找到每一个元胞数组中与输入的bin接近的所有值
        % 参数:
        % -- cellArray：元胞数组
        % -- ampArray: 元胞数组对应的能量数组
        % -- bin：要查找的bin
        % 结果:
        % -- result: 每一个元胞数组中与输入的bin接近的所有值
        % -- ampResult: 每一个元胞数组中与输入的bin接近的所有值对应的能量
        %% √
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

        %% 方法: 在元胞数组（n*m的矩阵）中找到满足syncword特征的第一个窗口位置
        % 参数:
        % -- BinArray：元胞数组
        % -- preamble_end_pos：满足sync word特征的
        % 结果: Preamble_end_pos
        %% √
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

        %% 利用 Preamble(Base-upchirp) 和 SFD(Base-downchirp) 相反偏移的性质，计算 CFO 和 winoffset
        function obj = getcfoWinoff(obj)
            % 计算主峰的 CFO (需要补零操作, 为了更好地评估峰值频率，可以通过用零填充原始信号来增加分析窗的长度。这种方法以更精确的频率分辨率自动对信号的傅里叶变换进行插值)
            % 对 Preamble 阶段的 FFT 峰值进行排序，得到前 filter 的峰
            zeropadding_size = obj.loraSet.factor;                   % 设置补零的数量，这里的 decim 表示，补上 decim-1 倍窗口的零，计算 FFT 时一共是 decim 倍的窗口（decim+1）, e.g. 16
            d_sf = obj.loraSet.sf;
            d_bw = obj.loraSet.bw;
            dine = obj.loraSet.dine;
            fft_x = obj.loraSet.fft_x;
            % Preamble_num = obj.loraSet.Preamble_length;
            downchirp = obj.idealDownchirp;
            upchirp = obj.idealUpchirp;
            % filter_num = obj.loraSet.filter_num;
            % leakage_width1 = obj.loraSet.leakage_width1;    % 0.0050
            % leakage_width2 = obj.loraSet.leakage_width2;    % 0.9950
            preambleSignalTmp = obj.preambleSignal;
            preambleEndPosTmp = obj.preambleEndPos;
            preambleBinTmp = obj.preambleBin;

            dine_zeropadding = dine * zeropadding_size * 2 ^ (10 - d_sf);   % e.g. 16384 * 16 * 2 ^ (10 - 10) = 262144
            fft_x_zeropadding = fft_x * zeropadding_size * 2 ^ (10 - d_sf);  % e.g. 1024 * 16 * 2 ^ (10 - 10) = 16384

            % 获取最后一个preamble窗口的若干个峰值，找到最接近preambleBin的峰
            samples = preambleSignalTmp((preambleEndPosTmp - 1) * dine + 1 : (preambleEndPosTmp) * dine);  % 倒数第二个preamble
            samples_fft = abs(fft(samples .* downchirp, dine_zeropadding, 2));
            samples_fft_merge = samples_fft(1 : fft_x_zeropadding) + samples_fft(dine_zeropadding - fft_x_zeropadding + 1 : dine_zeropadding);
            % 找出若干个峰值，峰值间间隔为leakWidth
            result = obj.findpeaksWithShift(samples_fft_merge, fft_x_zeropadding);
            binPos = result(2, :);
            % 找到与bin最接近的bin的索引
            %             closest_idx = obj.findClosetBin(binPos, preambleBin*zeropadding_size*2^(10-d_sf), fft_x_zeropadding);
            % [~, closest_idx] = min(abs(binPos - preambleBin*zeropadding_size*2^(10-d_sf))); % 找出最接近bin值的元素的索引
            % 找到与bin范围内接近的所有峰
            binArr = obj.findBinRange(binPos, preambleBinTmp * zeropadding_size * 2 ^ (10 - d_sf), zeropadding_size * 2, fft_x_zeropadding);
            [~, closest_idx] = min(abs(samples_fft_merge(binArr) - obj.preamblePeak));
            % upchirpBin = binPos(closest_idx); % 获得最接近preambleBin的bin作为用来对齐的upchirpBin
            upchirpBin = binArr(closest_idx);
            upchirpPeak = samples_fft_merge(upchirpBin); % 获取用于判断能量的标准值

            % 已知downchirp的位置，得到downchirp的bin（默认downchirp窗口内不存在downcrhip间的冲突）
            % TODO：可能需要考虑downcrhip发生冲突的问题
            SFD_samples = preambleSignalTmp((preambleEndPosTmp + 3) * dine + 1 : (preambleEndPosTmp + 4) * dine);  % 第一个和第二个downchirp之间
            SFD_samples_fft = abs(fft(SFD_samples .* upchirp, dine_zeropadding));
            samples_fft_merge = SFD_samples_fft(1 : fft_x_zeropadding) + SFD_samples_fft(dine_zeropadding - fft_x_zeropadding + 1 : dine_zeropadding);
            % 找出若干个峰值，峰值间间隔为leakWidth
            result = obj.findpeaksWithShift(samples_fft_merge, fft_x_zeropadding);
            if isempty(result)
                disp("找不到匹配的downchirp峰值");
                return;
            end
            [~, closest_idx] = min(abs(result(1, :) - upchirpPeak)); % 找出峰值能量最接近的元素的索引
            downchirpBin = result(2, closest_idx);

            % 计算CFO和窗口偏移量
            if upchirpBin + downchirpBin < fft_x_zeropadding * 0.5
                cfo_bin = upchirpBin + downchirpBin - 2;
                obj.cfo = -cfo_bin / 2 / fft_x_zeropadding * d_bw;
                obj.winOffset = (downchirpBin - upchirpBin) / 2 ^ (11 - d_sf);
            elseif upchirpBin + downchirpBin > fft_x_zeropadding * 1.5
                cfo_bin = upchirpBin + downchirpBin - fft_x_zeropadding*2 - 2;
                obj.cfo = -cfo_bin / 2 / fft_x_zeropadding * d_bw;
                obj.winOffset = (downchirpBin - upchirpBin) / 2 ^ (11 - d_sf);
            else
                cfo_bin = upchirpBin + downchirpBin - fft_x_zeropadding - 2;
                obj.cfo = -cfo_bin / 2 / fft_x_zeropadding * d_bw;
                obj.winOffset = (fft_x_zeropadding - (upchirpBin - downchirpBin)) / 2 ^ (11 - d_sf);
            end
        end


        %% 方法: 检测该信道该窗口内是否存在信号
        % 参数:
        % -- signals: 信号
        % 结果: isActive：布尔值，代表是否存在信号
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

        %%
        function binIndex = findBin(obj, binPos, findBin, fft_xTmp)
            arrLength = length(binPos);
            binPosRing = [binPos - fft_xTmp, binPos, binPos + fft_xTmp];
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

        %%
        function binPosRing = findBinRange(obj, binPos, findBin, range, fft_xTmp)
            binPosRing = [binPos - fft_xTmp, binPos, binPos + fft_xTmp];
            binPosRing = binPosRing(binPosRing >= findBin - range & binPosRing <= findBin + range);
            % [~, closestIndex] = min(abs(binPosRing - findBin));
            for index = 1:length(binPosRing)
                if binPosRing(index) <= 0
                    binPosRing(index) = binPosRing(index) + fft_xTmp;
                elseif binPosRing(index) > fft_xTmp
                    binPosRing(index) = binPosRing(index) - fft_xTmp;
                end
            end
        end

        %%
        function binIndex = findClosetBin(obj, binPos, findBin, fft_xTmp)
            arrLength = length(binPos);
            binPosRing = [binPos - fft_xTmp, binPos, binPos + fft_xTmp];
            [~, closestIndex] = min(abs(binPosRing - findBin));
            if closestIndex <= arrLength
                binIndex = closestIndex;
            elseif closestIndex > arrLength*2
                binIndex = closestIndex - 2*arrLength;
            else
                binIndex = closestIndex - arrLength;
            end
        end

%         function obj = getSFDPos(obj)
%             fft_xTmp = obj.loraSet.fft_x;
%             dine = obj.loraSet.dine;
%             preamble_len = obj.loraSet.Preamble_length;
%             binRecordSFD = zeros(1, preamble_len + 4);

% %             figure(1);
% %             FFT_plot(obj.preambleSignal(1:(preamble_len+4)*dine), obj.loraSet, obj.cfoDownchirp, 12);

%             for windows = 1 : preamble_len + 4
%                 signal = obj.preambleSignal((windows - 1) * dine + 1 : windows * dine);
%                 dechirp = signal .* obj.cfoDownchirp;
%                 dechirp_fft = abs(fft(dechirp, dine));
%                 dechirp_fft = dechirp_fft(1 : fft_xTmp) + dechirp_fft(dine - fft_xTmp + 1 : dine);
%                 [~, binRecordSFD(windows)] = max(dechirp_fft);
%             end
%             for windows = 3 : preamble_len + 4
%                 first = binRecordSFD(windows - 2);
%                 second = binRecordSFD(windows - 1);
%                 third = binRecordSFD(windows);
%                 if  (second - first >= 7 && second - first <= 9) && (third - second >= 7 && third - second <= 9)
%                     SFDPosTemp = windows;
%                     break;
%                 end
%             end
%             % 处理找不到最后一个 preamble 的情况
%             if exist('SFDPos', 'var') == 0
%                 SFDPosTemp = 10;
%             end
%             obj.SFDPos = SFDPosTemp;
%         end

        %% 方法: 获取 SFD 的位置
        %{
            因为信号窗口被调整了，所以需要重新寻找SFD的位置，保证正确率
            preambleSignal：调整窗口后已划分信道的信号
        %}
        function obj = getSFDPos(obj)
            fft_xTmp = obj.loraSet.fft_x;
            dine = obj.loraSet.dine;
            % leakWidth = obj.loraSet.leakage_width1;
            % preambleBin = obj.preambleBin;
            preambleEndPosTmp = obj.preambleEndPos;
            preambleSignalTmp = obj.preambleSignal;
            if preambleEndPosTmp < 2
                disp("preambleEndPos位置无效");
                return;
            end

            % 用于存储preambleLength+4数目窗口下，每个窗口的峰值
            candidate = cell(1, 5);
            candidatePeaks = cell(1, 5);

%             FFT_plot(preambleSignalTmp, obj.loraSet, obj.cfoDownchirp, 5);
            % 同样思路，找到这几个窗口内的前20峰值，然后找到与bin1最接近的峰，判断其规律是否与sync word相同
            for windows = 1 : 5
%                 signal = preambleSignalTmp((windows-1+preambleEndPosTmp-preamble_len)*dine+1 : (windows+preambleEndPosTmp-preamble_len)*dine);
                signal = preambleSignalTmp((windows - 1 + preambleEndPosTmp - 2) * dine + 1 : (windows + preambleEndPosTmp - 2) * dine);
                dechirp = signal .* obj.cfoDownchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1 : fft_xTmp) + dechirp_fft(dine - fft_xTmp + 1 : dine);
                result = obj.findpeaksWithShift(dechirp_fft, fft_xTmp);

                % 记录峰值
                candidate{windows} = result(2, :);
                candidatePeaks{windows} = result(1, :);
            end

            % 找到每个窗口中，最接近 preambleBin 的值
            [BinArray, AmpArray] = obj.findClosetSyncWordBin(candidate, candidatePeaks, 1);

            % 根据sync word的特征找到最后一个preamble
            SFDPosTmp = obj.findPosWithSyncWord(BinArray, AmpArray);

            % 处理找不到最后一个preamble，报错
            if exist('SFDPosTmp', 'var') == 0 || SFDPosTmp == 0
                disp("找不到SFD的位置");
                return;
            else
                obj.SFDPos = SFDPosTmp + preambleEndPosTmp - 2;
                % 获取SFD的能量，为后面找downchirp sync做准备
                peak = zeros(1, 2);
                obj.SFDPeakAmp = zeros(1, 2);
                for t = 0 : 1
                    signal = preambleSignalTmp((obj.SFDPos + t) * dine + 1 : (obj.SFDPos + t + 1 ) * dine);
                    dechirp = signal .* obj.cfoUpchirp;
                    dechirp_fft = abs(fft(dechirp, dine));
                    dechirp_fft = dechirp_fft(1 : fft_xTmp) + dechirp_fft(dine - fft_xTmp + 1 : dine);
                    % 冲突情况下，需要找到downchirp bin为1的峰的峰值
                    result = obj.findpeaksWithShift(dechirp_fft, fft_xTmp);
                    binPos = result(2, :);
                    closest_idx = obj.findClosetBin(binPos, 1, fft_xTmp);
                    if isempty(closest_idx)
                        disp("找不到SFD中满足要求的峰");
                        return;
                    end
                    peak(t+1) = result(1, closest_idx);
                    obj.SFDPeakAmp(t + 1) = peak(t + 1);
                end
                obj.peakStandard = mean(peak);
            end
        end
    end
end