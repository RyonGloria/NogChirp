% 将输入文件夹中的信号文件按信道读取
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 采样值文件读取路径和保存路径
readDir = 'E:\share\samples\';
% outDir = 'E:\tmp\';
% 读取配置和验证文件
[loraSet] = readLoraSet('sf7_BW250.json');

% 生成preamble对应信道的idealchirp
[downchirp, upchirp] = buildIdealchirp(loraSet, 0); % build idealchirp
% 读取文件夹下所有采样值文件
fileIn = dir(fullfile(readDir, '*.sigmf-data'));          
for fileCount = 1:1
    if mod(fileCount, 100) == 0           % 每循环100次，输出当前进度
        fprintf("The time is %d\n", fileCount);
    end
    % 从文件中读取信号
    [signal] = readSignalFile(readDir, fileIn(fileCount));
    % 将信号按照信道存放规则划分
    signal = reshape(signal, [], 4).';
%     signal = sum(signal, 1);
    figure(1);
    stft(signal(1, 1:20*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
    figure(2);
    stft(signal(2, 1:20*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
    figure(3);
    stft(signal(3, 1:20*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
end