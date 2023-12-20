function [] = writeBinChannel(loraSet, chirp1Bin, channelPredictArray, writeDir, PreambleChannel, channelNum, subchirpNum, downchirpSync, cfo)
    writePath = strcat(writeDir, 'sf', string(loraSet.sf), '\channel', string(channelNum), '\preamble', string(PreambleChannel), '\subchirpNum', string(subchirpNum), '\');
    if ~exist(writePath, 'dir')
        mkdir(writePath);
    end
    writeBinPath = strcat(writePath, 'bin_', 'downchirpsync', string(downchirpSync), '.txt');
    writeChannelPath = strcat(writePath, 'channel_', 'downchirpsync', string(downchirpSync), '.txt');
    writeCfoPath = strcat(writePath, 'cfo_', 'downchirpsync', string(downchirpSync), '.txt');
    writematrix(chirp1Bin, writeBinPath);
    writematrix(channelPredictArray, writeChannelPath);
    writematrix(cfo, writeCfoPath);
    fclose all;