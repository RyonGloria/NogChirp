classdef NogChirpDecoder < LoraDecoder
    properties
        record;
        payloadBin;
        binRecord;
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
        function obj = NogChirpDecoder(loraSet)
            obj@LoraDecoder(loraSet);
        end

        function obj = decode(obj, signals)
            obj.binRecord = [];
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

            obj = obj.singleDecode();
        end

        function obj = singleDecode(obj)
            fftX = obj.loraSet.fft_x;
            dine = obj.loraSet.dine;
            preambleSignalTemp = obj.preambleSignal;
            binArr = zeros(1, obj.loraSet.payloadNum);

            % dechirp做FFT，获得最大峰值即为信道矩阵的信息
            signal = preambleSignalTemp((obj.SFDPos + 2.25)*dine + 1 : end);
            for window = 1 : obj.loraSet.payloadNum
                signals = signal((window - 1) * dine + 1 : window * dine);
                dechirp = signals .* obj.cfoDownchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1 : fftX) + dechirp_fft(dine - fftX + 1 : dine);
                [~, binArr(window)] = max(dechirp_fft);
            end

            obj.binRecord = binArr;
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
            obj.PreambleBin = Preamble_bin;   % e.g. 823
        end

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