function [cfo, signalOut, packageFlag] = alignSignal(loraSet, signal, downchirp, upchirp, preambleChannel)
    packageFlag = true;
    [signalToProcess] = signalFrequencyShift(loraSet, signal, preambleChannel);
    [signalToProcess] = lowPassFilterFir(signalToProcess, loraSet);
    % find Preamble
    [PreambleStartPos] = detect_preamble(signalToProcess, loraSet, downchirp); 
    if PreambleStartPos == 999
        packageFlag = false;
    elseif PreambleStartPos ~= 1
        signal = circshift(signal, -(PreambleStartPos-1) * loraSet.dine);
    end
    % 调整cfo和timeoffset
    [cfo, windowsOffset] = get_cfo_winoff(signalToProcess, loraSet, downchirp, upchirp, loraSet.factor, true);
    signalOut = circshift(signal, -round(windowsOffset));