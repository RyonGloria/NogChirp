function [channel_jump_array] = getDownchirpSync(loraSet, signal, d_upchirp_cfo)
%     fft_x = loraSet.fft_x;
%     dine = loraSet.dine;
%     downchirpSyncSignal = signal(10*loraSet.dine+1:11*loraSet.dine);
%     dechirp = downchirpSyncSignal .* upchirp;
%     dechirp_fft = abs(fft(dechirp));
%     dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine);
%     [~, downchirpSync] = max(dechirp_fft);

    fftX = loraSet.fft_x;
    dine = loraSet.dine;
    % get downchirp_sync bin
%     downchirp_sync_signal = signal_detect(14.25*loraSet.dine+1:15.25*loraSet.dine);
    dechirp = signal .* d_upchirp_cfo;
    dechirp_fft = abs(fft(dechirp, dine));
    dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
    [tmp, downchirp_bin] = max(dechirp_fft);
    load(".\Result\random_record.mat");
    channel_jump_array = random_record(downchirp_bin, :);