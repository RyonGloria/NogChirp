function [signalOut] = signalFrequencyShift(loraSet, signal, carrirFre)
    Fs = loraSet.sample_rate;
    t = 0:1/Fs:1/Fs*(loraSet.dine-1);
    signalLength = length(signal);
    chirpNum = ceil(signalLength/loraSet.dine);
    m = repmat(exp(1i*2*pi*carrirFre*t), 1, chirpNum);
    signalOut = 2 .* signal .* m;