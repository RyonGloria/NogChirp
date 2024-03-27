classdef VarCutChirpDecoder < LoraDecoder
    properties
        payloadBin;
        BinRecord;  % 记录临时的bin值
        binRm;      % 记录目标数据包间隔内需要排除的其它 preamble 的 bin
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
        detectedPktAll;        % 记录解码信号的信息 (Channel, preambleBin, preambleEndPos, preamblePeak, SFDPos, CFO, winOffset, rebuildSignal)
        detectArray;           % 存放检测到的 bin
        detectArrayCount;      % 计数器
        detectArrayNum;        % 记录队列中存放有效数据的数目
    end

    methods
        %% 初始化方法
        function obj = VarCutChirpDecoder(loraSet, DebugUtil)
            obj@LoraDecoder(loraSet);
            % 初始化 debug 工具
            obj.DebugUtil = DebugUtil;
            obj.DebugUtil.info("", "初始化 VarCutChirpDecoder");
        end

        %% 方法: 清空类中的某些中间值
        function obj = clear(obj)
            % 清空类中的某些中间值
            obj.BinRecord = {};
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
            obj.errorFlag = false; % 每次进行setArgs相当于一次初始化，所以要对errorFlag置为false
            obj.errorMsg = "";
            obj.detectArray = zeros(1, obj.loraSet.fft_x);   % 存放信道检测到的峰值bin（bin 最多为 fft_x 种）
            obj.detectArrayCount = zeros(1, obj.loraSet.fft_x);  % 计数器
            obj.detectArrayNum = 0;  % 记录信道队列中存放有效数据的数目
        end

        %% 方法: 解码
        function obj = decode(obj, signals)
            obj = obj.clear();
            dineTmp = obj.loraSet.dine;
            signalsTmp = signals;
            obj.payloadBin = {};
            obj.detectedPktAll = {};

            % obj.channelFlag = true;
            %% 信道检测, 得到信道数以及每个信道中检测到的冲突数据包
            for Ch_i = 1 : 2
                % disp("Channel -" + num2str(Ch_i) + "-");
                if Ch_i == 2
                    signalsTmp = obj.rechangeSignalFreq(signals, - obj.loraSet.bw / 4);  % 中心频率对齐第二信道
                end

                windowsNum = fix(length(signalsTmp) / dineTmp); % 计算需要移动多少个窗口进行信号检测和解调
                % FFT_plot(obj.preambleSignal, obj.loraSet, obj.idealDownchirp, 23);
                for window_i = 1 : windowsNum % 开始扫描每一个窗口，发现其中的 premable
                    windowChirp = signalsTmp((window_i - 1) * dineTmp + 1 : window_i * dineTmp);
                    [obj, preambleBinTmp] = obj.detect(windowChirp);
                    if ~isempty(preambleBinTmp)
                        % obj.DebugUtil.debug("\t", "♦ preamble bin: [" + regexprep(num2str(preambleBinTmp), '\s+', ' ') + "]");
                        % disp("♦ preamble bin: [" + regexprep(num2str(preambleBinTmp), '\s+', ' ') + "]");
                        % 根据detect的结果进行解调
                        for num = 1 : length(preambleBinTmp)
                            obj = obj.clear();

                            obj.window = window_i;
                            obj.preambleBin = preambleBinTmp(num);
                            obj.DebugUtil.debug("\t", " preambleBin: " + obj.preambleBin);
                            obj.preambleSignal = signalsTmp;

                            %% 检测 preamble 结束位置
                            obj = obj.detectPreambleEndPosBehind();
                            if obj.errorFlag == true
                                continue;
                            end
                            obj.DebugUtil.debug("\t", " preambleEndPos: " + obj.preambleEndPos);
                            obj.DebugUtil.debug("\t", " preamblePeak: " + obj.preamblePeak);


                            %% 通过 preamble 和 SFD 的 bin 来计算 CFO 和 winoffset
                            obj = obj.getcfoWinoff();
                            if obj.errorFlag == true
                                continue;
                            end
                            obj.DebugUtil.debug("\t", " CFO: " + obj.cfo);
                            obj.DebugUtil.debug("\t", " winOffset: " + obj.winOffset);

                            %% 根据 winoffset 调整信号
                            obj.preambleSignal = circshift(obj.preambleSignal, -round(obj.winOffset));
                            obj = obj.rebuildIdealchirpCfo(0);  % 根据 cfo 重新生成带有 decfo 的 idealchirp，用于解调

                            %% 获取 SFD 的位置
                            obj = obj.getSFDPos();   % 获取 SFD 的位置
                            if obj.errorFlag == true
                                continue;
                            end
                            % disp("  • SFD Position: " + num2str(obj.SFDPos));
                            obj.DebugUtil.debug("\t", " SFD Real Position: " + obj.SFDPos + "\n");
                            % 信号所在信道, preamble Bin, preamble 结束位置, preamble 峰值, SFD 准确位置, CFO, winOffset, 修正 winOffset 后的信号
                            readytoAdd = {Ch_i, obj.preambleBin, obj.preambleEndPos, obj.preamblePeak, obj.SFDPos, obj.cfo, obj.winOffset, obj.preambleSignal};
                            obj.detectedPktAll = [obj.detectedPktAll; readytoAdd];
                        end
                    end
                end
            end

            %% 去除目标信号数据包中包含其它冲突信号的 preamble 的 bin 值 (因为冲突的 preamble 是 8 个连续的一样的 bin 值，在目标信号的数据包中产生的干扰显著，表现为占满整个解调窗口的错误 Bin，影响目标 chirp 的解调)
            numRows = size(obj.detectedPktAll, 1);
            obj.binRm = cell(1, numRows);   % 创建一个cell数组
            preamble_len = obj.loraSet.Preamble_length;
            for i = 1:numRows
                rowofBinRm = cell(1, obj.loraSet.payloadNum);  % 存放每一行中需要 remove 的 preamble bin
                targetPktStartPos = round(obj.detectedPktAll{i, 5} + 2.25);
                cfoTmp = obj.detectedPktAll{i, 6};          % 目标信号的 cfo
                winOffsetTmp = obj.detectedPktAll{i, 7};    % 目标信号的 winOffset
                for j = 1:numRows
                    if i == j
                        continue;
                    end
                    otherPreambleStartPos = round(obj.detectedPktAll{j, 5} - 10);
                    diff_val = otherPreambleStartPos - targetPktStartPos;
                    overlapStartPos = max(diff_val + 1, 1);                % 重叠的开始位置
                    overlapEndPos = min(diff_val + preamble_len + 1, 23);  % 重叠的结束位置
                    overlapIndices = overlapStartPos:overlapEndPos;        % 重叠的位置
                    binTmp = obj.detectedPktAll{j, 2};                     % 冲突信号的 preamble bin
                    otherBin = mod(round(binTmp + cfoTmp / obj.loraSet.bw * obj.loraSet.fft_x ...
                             + round(winOffsetTmp) / obj.loraSet.dine * obj.loraSet.fft_x ...
                             + 0.25 * obj.loraSet.fft_x), obj.loraSet.fft_x); % 冲突信号在目标信号解码窗口的 preamble bin
                    rowofBinRm(overlapIndices) = cellfun(@(x) [x otherBin], rowofBinRm(overlapIndices), 'UniformOutput', false);
                end
                obj.binRm{i} = rowofBinRm;
            end

            %% 根据上面检测到的信道和冲突数据包信息，进行解码
            for pkt_i = 1 : numRows
                obj = obj.clear();
                Ch_i = obj.detectedPktAll{pkt_i, 1};                 % 信道
                obj.preamblePeak = obj.detectedPktAll{pkt_i, 4};     % preamble 峰值
                obj.SFDPos = obj.detectedPktAll{pkt_i, 5};           % SFD 位置
                obj.cfo = obj.detectedPktAll{pkt_i, 6};              % cfo
                obj.winOffset = obj.detectedPktAll{pkt_i, 7};        % winOffset
                obj.preambleSignal = obj.detectedPktAll{pkt_i, 8};   % 修正 winOffset 后的信号
                obj = obj.rebuildIdealchirpCfo(0);  % 根据 cfo 重新生成带有 decfo 的 idealchirp，用于解调

                obj = obj.OverlapChDecoder(Ch_i, pkt_i);
                obj.payloadBin{end + 1} = obj.BinRecord;
            end
        end

        %% 方法: 重叠信道解码
        % 参数:
        % 结果: obj.payloadBin
        %% ✔
        function obj = OverlapChDecoder(obj, ChNum, pkt_i)
            dineTmp = obj.loraSet.dine;
            preambleSignalTemp = obj.preambleSignal;
            binRmTmp = obj.binRm{pkt_i};  % 需要去除的冲突信号的 preamble bin

            signal = preambleSignalTemp((obj.SFDPos + 2.25) * dineTmp + 1 : end);  % 仅包含数据包部分的目标信号

            for window_i = 1 : obj.loraSet.payloadNum
                windowChirp = signal((window_i - 1) * dineTmp + 1 : window_i * dineTmp);  % symbol

                windowChirp = obj.filterOutOtherCH(windowChirp, 150, 0, obj.loraSet.bw / 2);    % 只保留目标信号所在的信道（针对重叠信道问题）
                %% 对齐窗口能量方差
                obj = obj.decodeVarofPower(windowChirp, 1/16, 16, binRmTmp{window_i});    % 信号, 切分窗口大小, 步长, 错误 Bin

            end
        end

        %% 方法: 能量方差解码
        % 参数:
        % -- chirp: 窗口信号
        % 结果: obj.BinRecord
        %% √
        function obj = decodeVarofPower(obj, chirp, winSize, step, binRmTmp)
            fft_xTmp = obj.loraSet.fft_x;
            dineTmp  = obj.loraSet.dine;

            dechirp_fft = obj.decodeChirp(chirp);
            signalPeakInfo = obj.findpeaksWithShift(dechirp_fft, fft_xTmp);   % 找峰值
            signalPeakInfo = obj.powerExtraction(signalPeakInfo, 1);          % 滤掉低于能量阈值（preamble 1/2能量）的峰值
            for i = 1 : length(binRmTmp)
                signalPeakInfo(:, signalPeakInfo(2, :) == binRmTmp(i)) = [];  % 去除其它包的 preamble 值
            end
            % disp(['Candicate Bins: [', num2str(signalPeakInfo(2,:)), ']']);

            numPeaks = size(signalPeakInfo, 2);
            if numPeaks == 0
                obj.BinRecord = [obj.BinRecord NaN];
                return;
            elseif numPeaks == 1    % 如果只有一个峰值, 则直接返回
                obj.BinRecord = [obj.BinRecord signalPeakInfo(2, numPeaks)];
                return;
            end

            WinNum = (fft_xTmp - (fft_xTmp * winSize)) / step;
            ResInfo = cell(1, WinNum);
            for i = 0 : 1: WinNum
                S_t = obj.downchirpSW(winSize, i * step) .* chirp;
                dechirp_fft = abs(fft(S_t, dineTmp));      % fft()函数的结果是个复数, 取绝对值表示幅值
                dechirp_fft = dechirp_fft(1 : fft_xTmp) + dechirp_fft(dineTmp - fft_xTmp + 1 : dineTmp);
                ResInfo{i + 1} = [dechirp_fft; 1 : fft_xTmp];     % 记录每个滑动窗口解码结果（{peak, bin}, {peak, bin}, ...}）, size is 'fft_xTmp'
            end

            varianceSet = zeros(1, numPeaks);
            % for j = 1 : numPeaks
            %     peakInfo = zeros(1, WinNum + 1);     % length(ResInfo) = WinNum + 1
            %     for i = 1 : WinNum + 1
            %         peakInfo(i) = ResInfo{i}(1, signalPeakInfo(2, j));
            %     end
            %     varianceSet(j) = var(peakInfo);
            % end
            for j = 1 : numPeaks
                peakInfo = cellfun(@(x) x(1, signalPeakInfo(2, j)), ResInfo);
                varianceSet(j) = var(peakInfo);
            end
            BinWithVariance = [varianceSet; signalPeakInfo(2, :)];  % [variance; bin]

            [~, index] = min(BinWithVariance(1, : ));
            obj.BinRecord = [obj.BinRecord BinWithVariance(2, index)];
        end

        %% 方法: 生成指定频率偏移的 downchirp
        % 参数:
        % -- Freq: 频率
        % 结果: [DownchirpOut]
        %% ✔
        function [DownchirpOut] = cfoFreqShiftDownchirp(obj, Freq)
            cmx = 1 + 1 * 1i;
            pre_dir = 2 * pi;
            Ns = obj.loraSet.bw / obj.loraSet.fft_x;  % 1s 内的 symbol 数 (Ns = 1 / Ts = bw / 2^SF)
            T = -0.5 * obj.loraSet.bw * Ns;
            d_dt = 1 / obj.loraSet.sample_rate;       % 采样点间间隔的时间
            t = d_dt * (0 : 1 : obj.loraSet.dine - 1);
            f0 = obj.loraSet.bw / 2 + obj.cfo + Freq;

            DownchirpOut = cmx * (cos(pre_dir .* t .* (f0 + T * t)) + sin(pre_dir .* t .* (f0 + T * t)) * 1i);
        end

        %% 方法: de-chirp
        % 参数:
        % -- chirp: 信号
        % -- downchirpPick: 下行 chirp 类型选择
        % -- Freq: 频率
        % 结果: [signalOut]
        %% ✔
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
        %% ✔
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
        %% ✔
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
        % 结果: groupInfo
        %% ✔
        function [groupInfo] = powerExtraction(obj, signalInfo, winSize)
            peak = signalInfo(1, :);
            pos = signalInfo(2, :);
            PowerAmp = obj.preamblePeak;  % 信号能量幅值, preamble 的峰值
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
        % 结果: downchirpSW
        %% ×
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
            dineTmp = obj.loraSet.dine;

            r_t = offset / obj.loraSet.fft_x;
            seg_1 = d_dt * (0 : r_t * (dineTmp - 1));
            seg_2 = d_dt * (r_t * (dineTmp - 1) : (r_t + winSize) * (dineTmp - 1));
            seg_3 = d_dt * ((r_t + winSize) * (dineTmp- 1) : dineTmp - 1);
            Downchirp_1 = noneArr(seg_1);
            Downchirp_2 = cmx * (cos(pre_dir .* seg_2 .* (f0 + T * seg_2)) + sin(pre_dir .* seg_2 .* (f0 + T * seg_2)) * 1i);
            Downchirp_3 = noneArr(seg_3);
            downchirpSW = [Downchirp_1, Downchirp_2, Downchirp_3];
        end

        %% 方法: 检测 preamble 的可能值
        % 参数:
        % -- chirp: 窗口信号
        % 结果: preambleBinTmp
        %% √
        function [obj, preambleBinTmp] = detect(obj, chirp)
            preambleBinTmp = [];
            % preambleBinTmp 用于存放检测出来的lora包的preamble bin
            dineTmp = obj.loraSet.dine;
            fft_x = obj.loraSet.fft_x;

            if isActive(obj, chirp)  % 如果这个信道存在信号
                % 对这个信道的信号乘以downchirp，做FFT，通过findpeaks函数找到突出峰值对应的bin，来判断preamble特性
                chirpTmp = chirp;
                dechirp = chirpTmp .* obj.idealDownchirp;
                dechirp_fft = abs(fft(dechirp, dineTmp));
                dechirp_fft = dechirp_fft(1 : fft_x) + dechirp_fft(dineTmp - fft_x + 1 : dineTmp);
                % 找出若干个峰值
                result = obj.findpeaksWithShift(dechirp_fft, fft_x);
                % deWinpeakposInfo = obj.powerExtraction(result, 1);  % 滤掉低于能量阈值的峰值
                binPos = result(2, :);
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
                countThreshold = 5;   % 阈值(即认为出现了 5 次以上的峰值才认为是有效的 preamble bin)
                for pos_i = 1:detectArrayNumTmp
                    if obj.detectArrayCount(countTmp) == countThreshold
                        % 将挑选出的lora包的preamble bin记录
                        preambleBinTmp = [preambleBinTmp obj.detectArray(countTmp)];
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

        %% 方法: 检测当前信道下对应 premableBin 的最后一个 preamble 的位置（为后续cfo计算准备）
        % TODO: 因为涉及到多解码，建议把 preambleEndPos 信息作为输出传出
        % 参数:
        % -- obj.preambleBin: preamble 的 bin 值
        % -- obj.preambleSignal: 信号
        % 结果: obj.preambleEndPos
        %% √
        function obj  = detectPreambleEndPosBehind(obj)
            dineTmp = obj.loraSet.dine;
            fft_xTmp = obj.loraSet.fft_x;
            preamble_len = obj.loraSet.Preamble_length;
            windowTmp = obj.window;
            preambleSignalTmp = obj.preambleSignal;
            preambleBinTmp = obj.preambleBin;

            % 用于存储 preambleLength + 4 数目窗口下，每个窗口的峰值
            candidate = cell(1, preamble_len);
            cancidateAmp = cell(1, preamble_len);

            % 搜寻第 7 个 preamble 前后的 chirp
            for t = windowTmp : windowTmp + 7
                % 每一个窗口做FFT后，将获得的结果找出若干个峰值
                signal_tmp = preambleSignalTmp(t * dineTmp + 1 : (t + 1) * dineTmp);
                dechirp = signal_tmp .* obj.idealDownchirp;
                dechirp_fft = abs(fft(dechirp, dineTmp));
                dechirp_fft = dechirp_fft(1 : fft_xTmp) + dechirp_fft(dineTmp - fft_xTmp + 1 : dineTmp);
                % 找出若干个峰值
                result = obj.findpeaksWithShift(dechirp_fft, fft_xTmp);
                % 记录峰值
                candidate{t - windowTmp + 1} = result(2, :);
                cancidateAmp{t - windowTmp + 1} = result(1, :);
            end

            % 找到每个窗口中，最接近 preambleBin 的值
            [BinArray, AmpArray] = obj.findClosetSyncWordBin(candidate, cancidateAmp, preambleBinTmp);

            % 根据 sync word 的特征找到最后一个 preamble
            Preamble_end_pos = obj.findPosWithSyncWord(BinArray, AmpArray);

            % 处理找不到的情况
            if exist('Preamble_end_pos', 'var') == 0 || Preamble_end_pos == 0
                obj.DebugUtil.warning("\t", " 找不到最后一个 preamble 的位置\n");
                obj.errorFlag = true;
                obj.errorMsg = "找不到最后一个 preamble 的位置";
                return;
            else
                obj.preambleEndPos = windowTmp + Preamble_end_pos - 2;
                % 获取 preamble 的能量
                signal_tmp = preambleSignalTmp((obj.preambleEndPos-1) * dineTmp + 1 : (obj.preambleEndPos) * dineTmp);
                dechirp = signal_tmp .* obj.idealDownchirp;
                dechirp_fft = abs(fft(dechirp, dineTmp));
                dechirp_fft = dechirp_fft(1 : fft_xTmp) + dechirp_fft(dineTmp - fft_xTmp + 1 : dineTmp);
                % 找出若干个峰值，峰值间间隔为 leakWidth
                result = obj.findpeaksWithShift(dechirp_fft, fft_xTmp);
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
            fft_xTmp = obj.loraSet.fft_x;
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
                        condition(i, k) = condition(i, k) + fft_xTmp;
                    elseif condition(i, k) > fft_xTmp
                        condition(i, k) = condition(i, k) - fft_xTmp;
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
            fft_xTmp = obj.loraSet.fft_x;
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
                            secondFindArr(1, value-6) = secondFindArr(1, value-6) + fft_xTmp;
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
                                    firstFindArr(1, value-6) = firstFindArr(1, value-6) + fft_xTmp;
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

        %% 方法: 利用 Preamble(Base-upchirp) 和 SFD(Base-downchirp) 相反偏移的性质，计算 CFO 和 winoffset
        % 参数:
        % -- preambleSignal: 信号
        % -- preambleBin: Preamble 的 bin 值
        % -- preambleEndPos: Preamble 的结束位置
        % 结果:
        % -- obj.cfo: CFO (频域) 偏移
        % -- obj.winOffset: 窗口 (时域) 偏移
        %% ✔
        function obj = getcfoWinoff(obj)
            % 计算主峰的 CFO (需要补零操作, 为了更好地评估峰值频率，可以通过用零填充原始信号来增加分析窗的长度。这种方法以更精确的频率分辨率自动对信号的傅里叶变换进行插值)
            % 对 Preamble 阶段的 FFT 峰值进行排序，得到前 filter 的峰
            zeropadding_size = obj.loraSet.factor;                   % 设置补零的数量，这里的 decim 表示，补上 decim-1 倍窗口的零，计算 FFT 时一共是 decim 倍的窗口（decim+1）, e.g. 16
            d_sf = obj.loraSet.sf;
            d_bw = obj.loraSet.bw;
            dineTmp = obj.loraSet.dine;
            fft_xTmp = obj.loraSet.fft_x;
            % Preamble_num = obj.loraSet.Preamble_length;
            downchirp = obj.idealDownchirp;
            upchirp = obj.idealUpchirp;
            % filter_num = obj.loraSet.filter_num;
            % leakage_width1 = obj.loraSet.leakage_width1;    % 0.0050
            % leakage_width2 = obj.loraSet.leakage_width2;    % 0.9950
            preambleSignalTmp = obj.preambleSignal;
            preambleEndPosTmp = obj.preambleEndPos;
            preambleBinTmp = obj.preambleBin;

            dine_zeropadding = dineTmp * zeropadding_size * 2 ^ (10 - d_sf);   % e.g. 16384 * 16 * 2 ^ (10 - 10) = 262144
            fft_x_zeropadding = fft_xTmp * zeropadding_size * 2 ^ (10 - d_sf);  % e.g. 1024 * 16 * 2 ^ (10 - 10) = 16384

            % 获取最后一个 preamble 窗口的若干个峰值，找到最接近preambleBin的峰
            samples = preambleSignalTmp((preambleEndPosTmp - 1) * dineTmp + 1 : (preambleEndPosTmp) * dineTmp);  % 倒数第二个preamble
            samples_fft = abs(fft(samples .* downchirp, dine_zeropadding, 2));
            samples_fft_merge = samples_fft(1 : fft_x_zeropadding) + samples_fft(dine_zeropadding - fft_x_zeropadding + 1 : dine_zeropadding);
            % 找出若干个峰值，峰值间间隔为leakWidth
            result = obj.findpeaksWithShift(samples_fft_merge, fft_x_zeropadding);
            binPos = result(2, :);
            % 找到与bin最接近的bin的索引
            %             closest_idx = obj.findClosetBin(binPos, preambleBin*zeropadding_size*2^(10-d_sf), fft_x_zeropadding);
            % [~, closest_idx] = min(abs(binPos - preambleBin*zeropadding_size*2^(10-d_sf))); % 找出最接近bin值的元素的索引
            % 找到与 bin 范围内接近的所有峰
            binArr = obj.findBinRange(binPos, preambleBinTmp * zeropadding_size * 2 ^ (10 - d_sf), zeropadding_size * 2, fft_x_zeropadding);
            if isempty(binArr)
                obj.DebugUtil.warning("\t", " 找不到与 bin 范围内接近的所有峰\n");
                obj.errorFlag = true;
                obj.errorMsg = "找不到与 bin 范围内接近的所有峰";
                return;
            end
            [~, closest_idx] = min(abs(samples_fft_merge(binArr) - obj.preamblePeak));
            % upchirpBin = binPos(closest_idx); % 获得最接近preambleBin的bin作为用来对齐的upchirpBin
            upchirpBin = binArr(closest_idx);
            upchirpPeak = samples_fft_merge(upchirpBin); % 获取用于判断能量的标准值

            % 已知downchirp的位置，得到downchirp的bin（默认downchirp窗口内不存在downcrhip间的冲突）
            % TODO：可能需要考虑downcrhip发生冲突的问题
            SFD_samples = preambleSignalTmp((preambleEndPosTmp + 3) * dineTmp + 1 : (preambleEndPosTmp + 4) * dineTmp);  % 第一个和第二个downchirp之间
            SFD_samples_fft = abs(fft(SFD_samples .* upchirp, dine_zeropadding));
            samples_fft_merge = SFD_samples_fft(1 : fft_x_zeropadding) + SFD_samples_fft(dine_zeropadding - fft_x_zeropadding + 1 : dine_zeropadding);
            % 找出若干个峰值，峰值间间隔为leakWidth
            result = obj.findpeaksWithShift(samples_fft_merge, fft_x_zeropadding);
            if isempty(result)
                obj.DebugUtil.warning("\t", " 找不到匹配的 downchirp 峰值\n");
                obj.errorFlag = true;
                obj.errorMsg = "找不到匹配的 downchirp 峰值";
                return;
            end
            [~, closest_idx] = min(abs(result(1, :) - upchirpPeak)); % 找出峰值能量最接近的元素的索引
            downchirpBin = result(2, closest_idx);

            % 计算 CFO 和窗口偏移量
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
            % 判断 CFO 是否超出范围
            if abs(obj.cfo) / obj.loraSet.bw >= 0.05
                obj.DebugUtil.warning("\t", " CFO 超出范围 (" + obj.cfo + ")\n");
                obj.errorFlag = true;
                obj.errorMsg = "CFO 超出范围";
                return;
            end
        end


        %% 方法: 检测该信道该窗口内是否存在信号
        % 参数:
        % -- signals: 信号
        % 结果: isActive：布尔值，代表是否存在信号
        %% √
        function [isActive] = isActive(obj, signals)
            dineTmp = obj.loraSet.dine;
            fft_xTmp = obj.loraSet.fft_x;
            isActive = false;
            % 通过FFT的峰值判断是否存在信号
            dechirp = signals .* obj.idealDownchirp;
            dechirp_fft = abs(fft(dechirp, dineTmp));
            dechirp_fft = dechirp_fft(1 : fft_xTmp) + dechirp_fft(dineTmp - fft_xTmp + 1 : dineTmp);
            % 找到最大值
            [amp, ~] = max(dechirp_fft);
            % 最大值超过窗口均值的十倍
            if amp > mean(dechirp_fft) * 10
                isActive = true;
            end
        end

        %%
        function binIndex = findBin(~, binPos, findBin, fft_xTmp)
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
        function binPosRing = findBinRange(~, binPos, findBin, range, fft_xTmp)
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
            binPosRing = [binPos - fft_xTmp, binPos, binPos + fft_xTmp];
            [~, closestIndex] = min(abs(binPosRing - findBin));
            arrLength = length(binPos);
            binIndex = mod(closestIndex - 1, arrLength) + 1;
        end

        %% 方法: 获取 SFD 的位置
        % 说明: 因为信号窗口被调整了，所以需要重新寻找SFD的位置，保证正确率
        % preambleSignal：调整窗口后已划分信道的信号
        % 结果: obj.SFDPos
        %% √
        function obj = getSFDPos(obj)
            fft_xTmp = obj.loraSet.fft_x;
            dineTmp = obj.loraSet.dine;
            % leakWidth = obj.loraSet.leakage_width1;
            % preambleBin = obj.preambleBin;
            preambleEndPosTmp = obj.preambleEndPos;
            preambleSignalTmp = obj.preambleSignal;
            if preambleEndPosTmp < 2
                obj.DebugUtil.warning("\t", " preambleEndPos 位置无效\n");
                obj.errorFlag = true;
                obj.errorMsg = "preambleEndPos 位置无效";
                return;
            end

            % 用于存储 preambleLength + 4 数目窗口下，每个窗口的峰值
            candidate = cell(1, 5);
            candidatePeaks = cell(1, 5);

            % 同样思路，找到这几个窗口内的前20峰值，然后找到与 bin1 最接近的峰，判断其规律是否与 sync word 相同
            for windows = 1 : 5
                signal = preambleSignalTmp((windows - 1 + preambleEndPosTmp - 2) * dineTmp + 1 : (windows + preambleEndPosTmp - 2) * dineTmp);
                dechirp = signal .* obj.cfoDownchirp;
                dechirp_fft = abs(fft(dechirp, dineTmp));
                dechirp_fft = dechirp_fft(1 : fft_xTmp) + dechirp_fft(dineTmp - fft_xTmp + 1 : dineTmp);
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
                obj.DebugUtil.warning("\t", " 找不到 SFD 的位置\n");
                obj.errorFlag = true;
                obj.errorMsg = "找不到 SFD 的位置";
                return;
            else
                obj.SFDPos = SFDPosTmp + preambleEndPosTmp - 2;
                % 获取 SFD 的能量，为后面找 downchirp sync 做准备
                peak = zeros(1, 2);
                obj.SFDPeakAmp = zeros(1, 2);
                for t = 0 : 1
                    signal = preambleSignalTmp((obj.SFDPos + t) * dineTmp + 1 : (obj.SFDPos + t + 1 ) * dineTmp);
                    dechirp = signal .* obj.cfoUpchirp;
                    dechirp_fft = abs(fft(dechirp, dineTmp));
                    dechirp_fft = dechirp_fft(1 : fft_xTmp) + dechirp_fft(dineTmp - fft_xTmp + 1 : dineTmp);
                    % 冲突情况下，需要找到downchirp bin为1的峰的峰值
                    result = obj.findpeaksWithShift(dechirp_fft, fft_xTmp);
                    binPos = result(2, :);
                    closest_idx = obj.findClosetBin(binPos, 1, fft_xTmp);
                    if isempty(closest_idx)
                        obj.DebugUtil.warning("\t", " 找不到 SFD 中满足要求的峰\n");
                        obj.errorFlag = true;
                        obj.errorMsg = "找不到 SFD 中满足要求的峰";
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