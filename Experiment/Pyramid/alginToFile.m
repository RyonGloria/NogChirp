fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 采样值文件读取路径和保存路径
inDir = 'D:\CHchirp_20Nodes\Pyramid_SF10_BW125K\Node1\';
writeDir = 'D:\CHchirp_20Nodes\Pyramid_SF10_BW125K\Node1\tmp\';
% 读取配置和验证文件
[loraSet] = readLoraSet('sf10_BW125.json', 10, 125e3, 2e6);
bin_ref = load('E:\CHchirp\Code\Config\Pyramid\SF10.txt')';
payloadNum = 18;
% 生成preamble对应信道的idealchirp
[downchirp, upchirp] = buildIdealchirp(loraSet, 0); % build idealchirp

% 读取文件夹下所有采样值文件
fileIn = dir(fullfile(inDir, '*.sigmf-data'));     
dine = loraSet.dine;
fftX = loraSet.fft_x;
for fileCount = 1:length(fileIn)
    if mod(fileCount, 100) == 0           % 每循环100次，输出当前进度
        fprintf("The time is %d\n", fileCount);
    end
    [signal] = readSignalFile(inDir, fileIn(fileCount));
    % find Preamble
    [PreambleStartPos] = detect_preamble(signal, loraSet, downchirp); 
    if PreambleStartPos == 999
        continue;
    elseif PreambleStartPos ~= 1
        signal = circshift(signal, -(PreambleStartPos-1) * loraSet.dine);
    end
    % 调整cfo和timeoffset
    [cfo, windowsOffset] = get_cfo_winoff(signal, loraSet, downchirp, upchirp, loraSet.Preamble_length, loraSet.factor, true);
    signal = circshift(signal, -round(windowsOffset));
    [downchirpCfo, upchirpCfo] = rebuild_idealchirp_cfo(loraSet, cfo, 0);
    bin = zeros(1,payloadNum);
    for binCount = 1:payloadNum
        chirp = signal((11.25+binCount)*dine+1: (12.25+binCount)*dine);
        dechirp = chirp .* downchirpCfo;
        dechirp_fft = abs(fft(dechirp));
        dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
        [~, bin(binCount)] = max(dechirp_fft);
    end

%     disp(bin);

    if(sum((bin-bin_ref) == 0) == payloadNum)
        signalLength = length(signal)*2;
        signalProcessed = zeros(1, signalLength);
        signalReal = real(signal);   signalImag = imag(signal);
        signalProcessed(1:2:signalLength-1) = signalReal;
        signalProcessed(2:2:signalLength) = signalImag;
    
        writePath = strcat(writeDir, 'signal_', string(fileCount), '.sigmf-data');
        fid=fopen(writePath, 'wb');
        fwrite(fid, signalProcessed, 'float32');
        fclose all;
    end

    
end

toc;
fclose all;