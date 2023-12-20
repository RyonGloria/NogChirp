fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

readDir = 'E:\Pyramid_samples\Samples_base\';
% writeDir = 'E:\Pyramid_samples\conflict\';
[loraSet] = readLoraSet('sf7_BW125.json');
fileIn = dir(fullfile(readDir, '*.sigmf-data')); 
bin_ref = load('E:\Pyramid_samples\bin_ref.txt')';
[downchirp, upchirp] = buildIdealchirp(loraSet, 0); % build idealchirp
% samples = zeros(1, 3*(40+12.25)*loraSet.dine);
% len = 2*(40+12.25)*loraSet.dine;
dine = loraSet.dine;
fftX = loraSet.fft_x;
for count = 1:length(fileIn)
    signal = readSignalFile(readDir, fileIn(count));
    [cfo, windowsOffset] = get_cfo_winoff(signal, loraSet, downchirp, upchirp, loraSet.factor, false);
    signal = circshift(signal, -round(windowsOffset));
    [downchirpCfo, upchirpCfo] = rebuild_idealchirp_cfo(loraSet, cfo, 0);
    bin = zeros(1,40);
    for binCount = 1:40
        chirp = signal((11.25+binCount)*dine+1: (12.25+binCount)*dine);
        dechirp = chirp .* downchirpCfo;
        dechirp_fft = abs(fft(dechirp));
        dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
        [~, bin(binCount)] = max(dechirp_fft);
    end
    if(sum((bin-bin_ref) == 0) ~= 40)
        fprintf(fileIn(count).name);
        fprintf("\n");
    end
%     fprintf("num is %d\n",sum((bin-bin_ref) == 0));
%     fprintf("num is %d\n",bin);
end