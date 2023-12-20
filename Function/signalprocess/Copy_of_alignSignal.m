function [cfo, signalOut, packageFlag] = Copy_of_alignSignal(loraSet, signal, downchirp, upchirp, preambleChannel)
    packageFlag = true;
    [signalToProcess] = signalFrequencyShift(loraSet, signal, preambleChannel);
%     stft(signalToProcess(1:40*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
    [signalToProcess] = lowPassFilter(loraSet, signalToProcess);
%     stft(signalToProcess(1:40*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
%     stft(downchirp, loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
%     stft(upchirp, loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
    % find Preamble
    [PreambleStartPos] = detect_preamble(signalToProcess, loraSet, downchirp); 
    if PreambleStartPos == 999
        packageFlag = false;
    elseif PreambleStartPos ~= 1
        signal = circshift(signal, -(PreambleStartPos-1) * loraSet.dine);
    end
    % 调整cfo和timeoffset
    [cfo, windowsOffset] = get_cfo_winoff(signal, loraSet, downchirp, upchirp, loraSet.factor, true);
    signalOut = circshift(signal, -round(windowsOffset));