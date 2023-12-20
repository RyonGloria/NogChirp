classdef LoraDecoderTmp < CHchirpDecoder
    properties
        binRecord;
    end

    methods
        function obj = LoraDecoderTmp(loraSet)
            obj@CHchirpDecoder(loraSet);
        end

        function obj = decode(obj, signals)
            obj.binRecord = [];
            obj.preambleSignal = signals;

            % 检测preamble，确定存在preamble并且获得第一个preamble出现的窗口和preamble数目
            obj = obj.detectPreambleBinBehind();

            % 通过preamble和SFD的bin来计算CFO和winoffset
            obj = obj.getcfoWinoff();

            % 调整信号的winoffset
            obj.preambleSignal = circshift(obj.preambleSignal, -round(obj.winOffset));

            % 根据cfo重新生成带有decfo的idealchirp，用于解调
            obj = obj.rebuildIdealchirpCfo(0);

            obj = obj.getSFDPos();

            obj = obj.singleDecode();
        end

        function obj = singleDecode(obj)
            fftX = obj.loraSet.fft_x;
            dine = obj.loraSet.dine;
            preambleSignal = obj.preambleSignal;
            binArr = zeros(1, obj.loraSet.payloadNum);

            % dechirp做FFT，获得最大峰值即为信道矩阵的信息
            signal = preambleSignal((obj.SFDPos+2.25)*dine+1 : end);
            for window = 1:obj.loraSet.payloadNum
                signals = signal((window-1)*dine+1 : window*dine);
                dechirp = signals .* obj.cfoDownchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
                [~, binArr(window)] = max(dechirp_fft);
            end

            obj.binRecord = binArr;
        end
    end
end