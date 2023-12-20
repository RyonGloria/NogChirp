% 将输入文件夹中的信号文件读取分信道提取出来写入目的文件夹中
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 设置preamble信道可选项
preambleChannelChoice = [-281.250e3, -93.750e3, 93.750e3, 281.250e3];
% 设置实验参数
PreambleChannel = 4;
payloadStartOffset = 13.25;
subchirpNum = 16;  % subchirp数目
payloadNum = 40;
channelNum = 4;
channelChoiceIndexTable = [1,2];
channelChoiceNum = length(channelChoiceIndexTable);  % 实际选择的信道数目
deleteFlag = true;
% 读取配置和验证文件
[loraSet] = readLoraSet('sf7_BW125.json');

% 生成0频段的idealchirp
[downchirp, upchirp] = buildIdealchirp(loraSet, 0);

dine = loraSet.dine;
fftX = loraSet.fft_x;
downchirpSync = circshift(downchirp, dine/2);
% samples = downchirpSync;
samples = downchirpSync + 10*upchirp;
samples = awgn(samples, -10, 'measured');
subplot(3,1,1);
stft(samples, loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);

subplot(3,1,2);
dechirp = samples .* upchirp;
dechirp_fft = abs(fft(dechirp));
dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
plot(dechirp_fft);

subplot(3,1,3);
cor = xcorr(samples, downchirp, dine);
cor = abs(cor(dine+1 : 2*dine));
plot(cor);

fclose all;