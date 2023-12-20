[loraSet] = readLoraSet('sf10_BW125.json', 2e6);
[downchirp, ~] = buildIdealchirp(loraSet, 0); % build idealchirp
[~, upchirp_1] = buildIdealchirp(loraSet, loraSet.bw/2); % build idealchirp
[~, upchirp_2] = buildIdealchirp(loraSet, loraSet.bw/2 + loraSet.bw/1024*3.5); % build idealchirp
dine = loraSet.dine;
fft_x = loraSet.fft_x;
% stft(upchirp, loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
upchirp = upchirp_1 + upchirp_2;
dechirp = upchirp .* downchirp;
dechirp_fft = abs(fft(dechirp, dine));
dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine);
figure(1);
plot(dechirp_fft);

dechirp = upchirp(1:dine/4) .* downchirp(1:dine/4);
dechirp_fft = abs(fft(dechirp, dine));
dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine);
figure(2);
plot(dechirp_fft);