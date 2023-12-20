function [signal] = readSignalFile(fileDir, fileDec)
    % 读取数据
    filePath = strcat(fileDir, fileDec.name);
    fid=fopen(filePath,'rb');
    [signalUnprocessed]=fread(fid,'float32')';
    length = size(signalUnprocessed, 2);
    signal = signalUnprocessed(1:2:length-1) + signalUnprocessed(2:2:length)*1i;
    fclose(fid);