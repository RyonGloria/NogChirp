fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器
[loraSet] = readLoraSet('sf10_BW125.json');
[downchirp, upchirp] = buildIdealchirp(loraSet, 0);
% stft(upchirp, loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
realTmp = real(upchirp);
imagTmp = imag(upchirp);
upchirp = imagTmp + realTmp*1i;
stft(upchirp, loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
