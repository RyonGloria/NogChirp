% 将输入文件夹中的信号文件读取分信道提取出来写入目的文件夹中
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 采样值文件读取路径和保存路径
inDir = 'E:\share\';
% 读取配置和验证文件
[loraSet] = readLoraSet('sf10_BW125.json', 'helloworld.txt');
channelChoiceMtrix = load('E:\tmp.txt');

% 生成preamble对应信道的idealchirp
% [downchirp, upchirp] = buildIdealchirp(loraSet, -187.5e3); % build idealchirp
[downchirp, upchirp] = buildIdealchirp(loraSet, 0);
% 读取文件夹下所有采样值文件
fileIn = dir(fullfile(inDir, '*.sigmf-data'));     
txtBin = dir(fullfile(inDir, '*.txt')); 
for fileCount = 1:10
    if mod(fileCount, 100) == 0           % 每循环100次，输出当前进度
        fprintf("The time is %d\n", fileCount);
    end
    [signal] = readSignalFile(inDir, fileIn(fileCount));
    % adjust cfo and timeoffset
    [cfo, signal, packageFlag] = alignSignal(loraSet, signal, downchirp, upchirp);
    if packageFlag == false
        continue;
    end
%     figure(1);
%     stft(signal(1:40*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(256),'OverlapLength',128,'FFTLength',loraSet.fft_x);

    % 将信号根据信道滤波划分
    [singalOut] = divideChannel(loraSet, signal, false);
    % 获得downchirp bin
    [downchirpCfo, upchirpCfo] = rebuild_idealchirp_cfo(loraSet, cfo, 0);
    [downchirpSync] = getDownchirpSync(loraSet, singalOut(1,:), upchirpCfo);
    fprintf('The pkg %d bin is: ',fileCount);
    fprintf(' %d ',downchirpSync);
    downchirpSync = downchirpSync - 1;
    channelChoice = mod(channelChoiceMtrix(downchirpSync, :), 4);

    fft_x = loraSet.fft_x;
    dine = loraSet.dine;
    Verbin = load(strcat(inDir, txtBin(fileCount).name));
    for chirpCount = 0:3
        chirpIntegrated = [singalOut(channelChoice(chirpCount*4+1)+1, (13.25+chirpCount)*dine+1 : (13.25+chirpCount)*dine+dine/4), ...
        singalOut(channelChoice(chirpCount*4+2)+1, (13.25+chirpCount)*dine+dine/4+1 : (13.25+chirpCount)*dine+dine/2), ...
        singalOut(channelChoice(chirpCount*4+3)+1, (13.25+chirpCount)*dine+dine/2+1 : (13.25+chirpCount)*dine+dine*3/4), ...
        singalOut(channelChoice(chirpCount*4+4)+1, (13.25+chirpCount)*dine+dine*3/4+1 : (13.25+chirpCount)*dine+dine)];
        dechirp = downchirpCfo .* chirpIntegrated;
        dechirp_fft = abs(fft(dechirp));
        dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine);
        [~, chirp1_bin] = max(dechirp_fft);
        fprintf(' %d ',chirp1_bin - Verbin(chirpCount+1));
    end
    figure(1);
    subplot(9,1,fileCount);
    plot(dechirp_fft);
    fprintf('\n');
    
    tmp = singalOut(1, 1:dine);
    dechirp = downchirpCfo .* tmp;
    dechirp_fft = abs(fft(dechirp));
    dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine);
    figure(2);
    subplot(9,1,fileCount);
    plot(dechirp_fft);
%     channel_choice = [0,1,1,3];
%     fft_x = loraSet.fft_x;
%     dine = loraSet.dine;
%     chirpIntegrated = [singalOut(channel_choice(1)+1, 13.25*dine+1 : 13.25*dine+dine/4), ...
%                         singalOut(channel_choice(2)+1, 13.25*dine+dine/4+1 : 13.25*dine+dine/2), ...
%                         singalOut(channel_choice(3)+1, 13.25*dine+dine/2+1 : 13.25*dine+dine*3/4), ...
%                         singalOut(channel_choice(4)+1, 13.25*dine+dine*3/4+1 : 13.25*dine+dine)];
%     
%     figure(3);
%     stft(chirpIntegrated, loraSet.sample_rate, 'Window',rectwin(256),'OverlapLength',128,'FFTLength',loraSet.fft_x);
%     figure(4);
%     stft(downchirp0, loraSet.sample_rate, 'Window',rectwin(256),'OverlapLength',128,'FFTLength',loraSet.fft_x);
%     
%     dechirp = downchirp0 .* chirpIntegrated;
%     dechirp_fft = abs(fft(dechirp));
%     dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine);
%     figure(5);
%     plot(dechirp_fft);
%     [~, chirp1_bin] = max(dechirp_fft);
%     disp(chirp1_bin);
end
fclose all;