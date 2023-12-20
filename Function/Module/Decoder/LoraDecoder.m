classdef LoraDecoder
    properties
        loraSet;
        idealUpchirp;
        idealDownchirp;
        cfoUpchirp;
        cfoDownchirp;
        cfo;
        winOffset;
    end
    
    methods
        function obj = LoraDecoder(loraSet)
            obj.loraSet = loraSet;
            obj = obj.buildIdealchirp(0);
        end
        
        function obj = buildIdealchirp(obj, f_temp)
            cmx = 1+1*1i;
            pre_dir = 2*pi;
            f0 = f_temp+obj.loraSet.bw/2;                           % 设置理想upchirp和downchirp的初始频率
            f1 = -f_temp+obj.loraSet.bw/2;
            d_symbols_per_second = obj.loraSet.bw / obj.loraSet.fft_x;
            T = -0.5 * obj.loraSet.bw * d_symbols_per_second;
            d_samples_per_second = obj.loraSet.sample_rate;        % sdr-rtl的采样率
            d_dt = 1/d_samples_per_second;         % 采样点间间隔的时间
            t = d_dt*(0:1:obj.loraSet.dine-1);
            % 计算理想downchirp和upchirp存入d_downchirp和d_upchirp数组中（复数形式）
            obj.idealDownchirp = cmx * (cos(pre_dir .* t .* (f0 + T * t)) + sin(pre_dir .* t .* (f0 + T * t))*1i);
            obj.idealUpchirp = cmx * (cos(pre_dir .* t .* (f1 + T * t) * -1) + sin(pre_dir .* t .* (f1 + T * t) * -1)*1i);
        end
        
        function obj = rebuildIdealchirpCfo(obj, f_temp)
            cmx = 1+1*1i;
            pre_dir = 2*pi;
            d_symbols_per_second = obj.loraSet.bw / obj.loraSet.fft_x;
            T = -0.5 * obj.loraSet.bw * d_symbols_per_second;
            d_samples_per_second = obj.loraSet.sample_rate;       % sdr-rtl的采样率
            d_dt = 1/d_samples_per_second;         % 采样点间间隔的时间
            t = d_dt*(0:1:obj.loraSet.dine-1);
            f0 = f_temp+obj.loraSet.bw/2+obj.cfo;                           % 设置理想upchirp和downchirp的初始频率
            f1 = -f_temp+obj.loraSet.bw/2-obj.cfo;
            
            % 计算理想downchirp和upchirp存入d_downchirp和d_upchirp数组中（复数形式）
            obj.cfoDownchirp = cmx * (cos(pre_dir .* t .* (f0 + T * t)) + sin(pre_dir .* t .* (f0 + T * t))*1i);
            obj.cfoUpchirp = cmx * (cos(pre_dir .* t .* (f1 + T * t) * -1) + sin(pre_dir .* t .* (f1 + T * t) * -1)*1i);
        end
        
        function obj = getcfoWinoff(obj, signals)
            % 计算主峰的CFO(需要补零操作)
            % 对Preamble阶段的FFT峰值进行排序，得到前filter的峰
            zeropadding_size = obj.loraSet.factor;                   % 设置补零的数量，这里的decim表示，补上decim-1倍窗口的零，计算FFT时一共是decim倍的窗口（decim+1）
            d_sf = obj.loraSet.sf;
            d_bw = obj.loraSet.bw;
            dine = obj.loraSet.dine;
            fft_x = obj.loraSet.fft_x;
            Preamble_num = obj.loraSet.Preamble_length;
            downchirp = obj.idealDownchirp;
            upchirp = obj.idealUpchirp;
            filter_num = obj.loraSet.filter_num;
            leakage_width1 = obj.loraSet.leakage_width1;
            leakage_width2 = obj.loraSet.leakage_width2;
            
            dine_zeropadding = dine*zeropadding_size*2^(10-d_sf);
            fft_x_zeropadding = fft_x*zeropadding_size*2^(10-d_sf);
            
            samples = reshape(signals(1:Preamble_num*dine),[dine,Preamble_num]).';  % dine 行 Preamble_num 列的矩阵
            samples_fft = abs(fft(samples .* downchirp, dine_zeropadding, 2));      % FFT 是沿着列方向进行的
            samples_fft_merge = samples_fft(:, 1:fft_x_zeropadding) + samples_fft(:,dine_zeropadding-fft_x_zeropadding+1:dine_zeropadding);
            [peak,pos] = sort(samples_fft_merge(1:Preamble_num,:),2,'descend');         % 对FFT进行排序
            Peak_pos = zeros(size(pos,1),filter_num);
            Peak_amp = zeros(size(peak,1),filter_num);
            Peak_pos(:,1) = pos(:,1);
            Peak_amp(:,1) = peak(:,1);
            for row = 1:size(pos,1)
                temp_array = ones(1,size(pos,2));
                for list = 1:filter_num
                    temp_array = temp_array & (abs(Peak_pos(row,list) - pos(row,:)) > fft_x_zeropadding*leakage_width1 & abs(Peak_pos(row,list) - pos(row,:)) < fft_x_zeropadding*leakage_width2);
                    temp_num = find(temp_array==1,1,'first');
                    Peak_pos(row,list+1) = pos(row,temp_num);
                    Peak_amp(row,list+1) = peak(row,temp_num);
                end
            end
            
            %寻找与第一个窗口的峰（默认第一个窗口只有包1的峰）相近的峰，得到与其相近且重复次数最多的bin，记作Preamble的bin
            if Peak_pos(1) == Peak_pos(2)
                upchirp_ref = Peak_pos(1);
            else
                upchirp_ref = Peak_pos(2);
            end
            upchirp_index = abs(Peak_pos-upchirp_ref) < fft_x_zeropadding*leakage_width1 | abs(Peak_pos-upchirp_ref) > fft_x_zeropadding*leakage_width2;
            upchirp_bin = (Peak_pos(upchirp_index));
            upchirp_peak = mode(upchirp_bin);
            
            % 已知downchirp的位置，得到downchirp的bin
            SFD_samples = signals((Preamble_num+2)*dine+1:(Preamble_num+3)*dine);
            SFD_samples_fft = abs(fft(SFD_samples .* upchirp, dine_zeropadding));
            samples_fft_merge = SFD_samples_fft(1:fft_x_zeropadding) + SFD_samples_fft(dine_zeropadding-fft_x_zeropadding+1:dine_zeropadding);
            [~, downchirp_peak] = max(samples_fft_merge);
            
            % 计算CFO和窗口偏移量
            if upchirp_peak + downchirp_peak < fft_x_zeropadding*0.5
                cfo_bin = upchirp_peak + downchirp_peak - 2;
                obj.cfo = -cfo_bin/2/fft_x_zeropadding * d_bw;
                obj.winOffset = (downchirp_peak - upchirp_peak) / 2^(11-d_sf);
            elseif upchirp_peak + downchirp_peak > fft_x_zeropadding*1.5
                cfo_bin = upchirp_peak + downchirp_peak - fft_x_zeropadding*2 - 2;
                obj.cfo = -cfo_bin/2/fft_x_zeropadding * d_bw;
                obj.winOffset = (downchirp_peak - upchirp_peak) / 2^(11-d_sf);
            else
                cfo_bin = upchirp_peak + downchirp_peak - fft_x_zeropadding - 2;
                obj.cfo = -cfo_bin/2/fft_x_zeropadding * d_bw;
                obj.winOffset = (fft_x_zeropadding - (upchirp_peak - downchirp_peak)) / 2^(11-d_sf);
            end
        end

        function [signalOut] = signalFrequencyShift(obj, signal, carrirFre)
            Fs = obj.loraSet.sample_rate;
            t = 0:1/Fs:1/Fs*(obj.loraSet.dine-1);
            signalLength = length(signal);
            chirpNum = ceil(signalLength/obj.loraSet.dine);
            m = repmat(exp(1i*2*pi*carrirFre*t), 1, chirpNum);
            signalOut = 2 .* signal .* m;
        end

        function [singalOut] = lowPassFilterFir(obj, signal)
            b = fir1(30, obj.loraSet.pass_arg, "low");
            singalOut = filter(b, 1, signal);
        end

        function [sortedGroups] = findpeaksWithShift(obj, fftResult, fft_x)
            % fft_x = obj.loraSet.fft_x;
            if obj.loraSet.sf == 10
                prominencesThreshold = 0.7;
            elseif obj.loraSet.sf == 9
%                 prominencesThreshold = 0.85;
                prominencesThreshold = 0.6;
            end
            leakWidth = obj.loraSet.leakage_width1;
%             [peak1, binPos1, ~, prominences1] = findpeaks(fftResult, "MinPeakDistance", fft_x*leakWidth, "SortStr", "descend");
            [peak1, binPos1, ~, prominences1] = findpeaks(fftResult, "SortStr", "descend");
            indecis = find((prominences1./peak1) >= prominencesThreshold);
            peak1 = peak1(indecis);
            binPos1 = binPos1(indecis);
            fftResultShift = circshift(fftResult, [0, fft_x/2]);  % 循环位移一半的窗口
%             [peak2, binPos2, ~, prominences2] = findpeaks(fftResultShift, "MinPeakDistance", fft_x*leakWidth, "SortStr", "descend");
            [peak2, binPos2, ~, prominences2] = findpeaks(fftResultShift, "SortStr", "descend");
            indecis = find((prominences2./peak2) >= prominencesThreshold);
            peak2 = peak2(indecis);
            binPos2 = binPos2(indecis);
            % 先对binPos2进行bin的矫正
            for i = 1:length(binPos2)
                binPos2(i) = binPos2(i) - fft_x/2;
                if binPos2(i) <= 0
                    binPos2(i) = binPos2(i) + fft_x;
                end
            end
            group1 = [peak1; binPos1];
            group2 = [peak2; binPos2];
            % 合并两组数组
            combinedGroups = [group1, group2];

            % 去重每组的第二个数组
            % 使用第二行数据进行去重
            [~, uniqueIndices, ~] = unique(combinedGroups(2, :), 'stable');

            % 提取去重后的数组
            uniqueGroups = combinedGroups(:, uniqueIndices);

            % 按照每组的第一个数组的大小进行排序
            sortedGroups = uniqueGroups;
            % 使用sortrows函数按第一行数据从大到小排序
            sortedGroups = sortrows(sortedGroups.', -1).';
        end
        
    end
end