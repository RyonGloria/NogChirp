function [] = writeSignalFile(loraSet, signal, writeDir, PreambleChannel, channelNum, subchirpNum, downchirpSync)
    signalLength = length(signal)*2;
    signalProcessed = zeros(1, signalLength);
    signalReal = real(signal);   signalImag = imag(signal);
    signalProcessed(1:2:signalLength-1) = signalReal;
    signalProcessed(2:2:signalLength) = signalImag;

    writePath = strcat(writeDir, 'sf', string(loraSet.sf), '\channel', string(channelNum), '\preamble', string(PreambleChannel), '\subchirpNum', string(subchirpNum), '\');
    if ~exist(writePath, 'dir')
        mkdir(writePath);
    end
    writePath = strcat(writePath, 'signal_', 'downchirpsync', string(downchirpSync), '.sigmf-data');
    fid=fopen(writePath, 'wb');
    fwrite(fid, signalProcessed, 'float32');
    fclose all;