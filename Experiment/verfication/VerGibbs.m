%% 
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 读取配置和验证文件
[loraSet] = readLoraSet('sf7_BW125.json');

% 生成0频段的idealchirp
dine = loraSet.dine;
fftX = loraSet.fft_x;
[downchirp, upchirp] = buildIdealchirp(loraSet, 0);

% samples = circshift(upchirp, 0);
% dechirp = samples .* downchirp;
% dechirp_fft = abs(fft(dechirp, dine));
% dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
% subplot(3,1,1);
% plot(dechirp_fft);
% 
% samples = circshift(upchirp, 16);
% dechirp = samples .* downchirp;
% dechirp_fft = abs(fft(dechirp));
% dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
% subplot(3,1,2);
% plot(dechirp_fft);
% 
% samples = circshift(upchirp, 0);
% dechirp = samples(1:dine/4) .* downchirp(1:dine/4);
% dechirp_fft = abs(fft(dechirp,dine));
% dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
% subplot(3,1,3);
% plot(dechirp_fft);

loraSet.channel = 0;
loraSet.channel_choice = 0;
loraSet.channel_choice_num = 0;
[d_upchirp, ~] = build_idealchirp_tmp(loraSet, 0);
[~, d_downchirp] = build_idealchirp_tmp(loraSet, 0);
dechirp = d_upchirp .* d_downchirp;
dechirp_fft = abs(fft(dechirp, dine));
dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
subplot(3,1,1);
plot(dechirp_fft);

[d_upchirp, ~] = build_idealchirp_tmp(loraSet, 64);
[~, d_downchirp] = build_idealchirp_tmp(loraSet, 0);
dechirp = d_upchirp .* d_downchirp;
dechirp_fft = abs(fft(dechirp, dine));
dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
subplot(3,1,2);
plot(dechirp_fft);

[d_upchirp, ~] = build_idealchirp_tmp(loraSet, 0);
[~, d_downchirp] = build_idealchirp_tmp(loraSet, 0);
dechirp = d_upchirp(1:dine/4) .* d_downchirp(1:dine/4);
dechirp_fft = abs(fft(dechirp, dine));
dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
subplot(3,1,3);
plot(dechirp_fft);