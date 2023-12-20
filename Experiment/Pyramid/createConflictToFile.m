% 画出所选文件的STFT，并给出bin值
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

readDir = 'E:\Pyramid_samples\Samples_base\';
writeDir = 'E:\Pyramid_samples\conflict\';
[loraSet] = readLoraSet('sf7_BW125.json');
fileIn = dir(fullfile(readDir, '*.sigmf-data')); 

% 根据参数生成发包时间矩阵
nodeNum = 10;
pkgRate = 4;
pktNum = 5;
pktLength = loraSet.dine * (12.25+40);
fs = 1e6;
[T] = gen_pkt_time_pyramid(nodeNum, pkgRate, pktNum, pktLength, fs);
load('C:\Users\ZKevin\Desktop\CHchirp\experiment\Pyramid\exp0\node_info.mat');

lengthAll = max(T,[],"all");
samples = zeros(1, ceil(lengthAll) + (40+12.25)*loraSet.dine);
T = T(:, 2:end);
T = T(:);
for count = 1:length(T)
    fail = node_info(ceil(1000*rand()), end);
    id = ceil(length(fileIn)*rand());
    signal = readSignalFile(readDir, fileIn(id));
    signal = sqrt(10^(fail/10)).*signal;
    samplesTmp = samples(ceil(T(count))+1 : ceil(T(count))+(40+12.25)*loraSet.dine);
    samplesTmp = samplesTmp + signal(1:(40+12.25)*loraSet.dine);
    samples(ceil(T(count))+1 : ceil(T(count))+(40+12.25)*loraSet.dine) = samplesTmp;
end
% 加入噪声
samples = awgn(samples, 20, 'measured');

% 写入文件中
signalLength = length(samples)*2;
signalProcessed = zeros(1, signalLength);
signalReal = real(samples);   signalImag = imag(samples);
signalProcessed(1:2:signalLength-1) = signalReal;
signalProcessed(2:2:signalLength) = signalImag;

writePath = strcat(writeDir, 'conflict.sigmf-data');
fid=fopen(writePath, 'wb');
fwrite(fid, signalProcessed, 'float32');


toc;
fclose all;