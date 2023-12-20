% 画出所选文件的STFT，并给出bin值
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 采样值文件读取路径和保存路径
inDir = 'E:\share\collision\';
% 读取配置和验证文件
[loraSet] = readLoraSet('sf10_BW125.json', 2e6);
bw = loraSet.bw;  % preamble占用的带宽
gbw = bw/2;     % 过渡带带宽
subchirpNum = 4;    % subchirp的数目
% 每个信道对应的中心频率（目前的信号已经进过频移至0频，其中心频率为采样时设置的中心频率）
preambleChannelChoice = [-(bw+gbw)/2-(bw+gbw), -(bw+gbw)/2, (bw+gbw)/2, (bw+gbw)/2+(bw+gbw)];
% 生成中心频率为0的idealchirp
[downchirp, upchirp] = buildIdealchirp(loraSet, 0); % build idealchirp
% 信号bin的groundtruth
true_bin = [0, 56, 112, 168, 224, 280, 336, 392, 448, 504, 560, 616, 672, 728, 784, 840] + 1;
% 读取文件夹下所有采样值文件
fileIn = dir(fullfile(inDir, '*.sigmf-data'));
% 检测队列
detect_array = zeros(4,40);   % 存放4个信道检测到的bin
detect_array_count = zeros(4,40);  % 计数器
detect_array_number = zeros(1,4);  % 记录四个信道队列中存放有效数据的数目
% preamble队列，存放在检测队列中连续出现了7次相同bin的峰值
preamble_array_bin = zeros(4,40);   % 存放4个信道确认preamble的bin
preamble_array_amp = zeros(4,40); 
preamble_array_count = zeros(4,40); % 计数器
preamble_array_number = zeros(1,4); % 同上，记录数目
% channel队列，在preamble队列中数4个窗口后计算downchirp_bin，
% cfo，winoff和跳信道策略矩阵对应的bin后放入到channel队列中
channel_array_cfo = zeros(4,40);  % 记录四个信道的cfo
channel_array_winoff = zeros(4,40); % 记录winoff
channel_array_amp = zeros(4, 40);
channel_array_count = zeros(4,40); % 计数器
channel_array_num = zeros(1,4); % 同上
% demodulate队列，在channel队列中一直等到第一个subchirp窗口进行重新对齐后，
% 放入demodulate队列中，在后续窗口中连续解调解码（不再区分信道）
demodulate_array_bin = zeros(160, 200);  % 存放解调的bin值
demodulate_array_channel = zeros(160, 200);  % 存放每一个数据包的跳信道策略
demodulate_array_cfo = zeros(1, 160);  % cfo
demodulate_array_winoff = zeros(1, 160);  % winoff
demodulate_array_count = zeros(1, 160);   % 计数器
demodulate_array_timeoff = zeros(1, 160); % timeoff
demodulate_array_num = 0;  % 同上
% demodulated队列，将demodulate队列中完成所有解码的数据放入，保存结果
demodulated_array = zeros(1600, 200);   % 结果数据
demodulated_array_num = 0; % 同上
for fileCount = 1:length(fileIn)
    % 从文件中读取信号流
    [signal] = readSignalFile(inDir, fileIn(fileCount));
    % 模拟gnuradio中的信号流
%     for windows = 1:floor(length(signal)/loraSet.dine)-1
    for windows = 1:100
        % 对每一个输入信号都进行correlate
        cor_value = detect_preamble_autocorr(signal((windows-1)*loraSet.dine+1 : (windows+1)*loraSet.dine), loraSet);
        % 存在preamble
        dine = loraSet.dine;
        fft_x = loraSet.fft_x;
        % 将信号划分成四个信道
        [signalOut] = divideChannel(loraSet, signal((windows-1)*loraSet.dine+1 : (windows+2)*loraSet.dine), preambleChannelChoice, false);
        % detect队列处理，后续测试发现，这个cor_value的方法不需要
        if(cor_value) > 0.01
            % 确定哪一个信道存在信号
            [active_channel] = get_active_channel_multi(signalOut, loraSet, downchirp);
            % 遍历每一个存在信号的信道，获得过滤的峰值，来确定preamble
            for channel_num = 1:size(active_channel,2)
                channel = active_channel(channel_num);
                signal_tmp = signalOut(channel, 1:dine);
                dechirp = signal_tmp .* downchirp;
                dechirp_fft = abs(fft(dechirp, dine));
                dechirp_fft = dechirp_fft(1:fft_x) + dechirp_fft(dine-fft_x+1:dine);
                % 对该窗口进行FFT，获得过滤取得前N个峰值，用作preamble的判定
                [pos_array, peak_array] = fft_filterN_peak(loraSet, dechirp_fft, 40);
                % 对峰值进行能量筛选（大于窗口均值的10倍）
                pos_pass = pos_array(peak_array > mean(dechirp_fft) * 10);
                % 对每一个筛选过的峰值进行记录
                for pos_count = 1:size(pos_pass, 2)
                    % 如果检测队列为空或者检测队列中没有该元素，直接加入
                    is_member = ismember(detect_array(channel, 1:detect_array_number(channel)), pos_pass(pos_count));
                    if detect_array_number(channel) == 0 || sum(is_member) == 0
                        % 检测队列的总数量加1
                        detect_array_number(channel) = detect_array_number(channel) + 1;
                        % 将峰对应的位置加入到检测队列中
                        detect_array(channel, detect_array_number(channel)) = pos_pass(pos_count);
                        % 并使其对应的count置1
                        detect_array_count(channel, detect_array_number(channel)) = 1;
                    else % 如果检测队列有该元素，对其记录的count加1
                        % 找到对应的检测队列中的位置
                        find_pos = find(detect_array(channel, 1:detect_array_number(channel)) == pos_pass(pos_count));
                        % 使其对应的count加1
                        detect_array_count(channel, find_pos) = detect_array_count(channel, find_pos) + 1;
                    end
                end
                % 将在检测队列中未能出现在此次preamble的峰值给剔除掉
                % 对检测队列中的每一个峰值进行遍历，找到没有在此次过滤出的峰出现的
                for pos_count = 1:detect_array_number(channel)
                    % 发现检测队里中的峰值在此次中未出现
                    is_member = ismember(pos_pass, detect_array(channel, pos_count));
                    if sum(is_member) == 0
                        % 在检测队列中进行删除
                        detect_array(channel, pos_count:end) = [detect_array(channel, pos_count+1:end), 0];
                        detect_array_count(channel, pos_count:end) = [detect_array_count(channel, pos_count+1:end), 0];
                        detect_array_number(channel) = detect_array_number(channel) - 1;
                    end
                end
                % 将检测队列中count数目等于7的元素放入检测成功的队列中
                for pos_count = 1:detect_array_number(channel)
                    if detect_array_count(channel, pos_count) == 7
                        % 重新进行补零后的fft，方便后续进行cfo和winoff的计算
                        [samples_fft_merge] = get_zeropadding_fft(signalOut(channel, 1:dine), loraSet, downchirp);
                        % 获得去掉旁瓣影响的前N个峰
                        [pos_array_padding, peak_array_padding] = fft_filterN_peak(loraSet, samples_fft_merge, 40);
                        % 找到原来没补零前的峰的位置
                        find_pos = find(detect_array(channel, pos_count) == pos_array);
                        % 找到补零fft后对应的峰值位置
                        pos = pos_array_padding(find_pos);
                        % 获取补零fft后对应的峰值强度
                        peak = peak_array_padding(find_pos);
                        % preamble队列对应信道元素数目加1
                        preamble_array_number(channel) = preamble_array_number(channel) + 1;
                        preamble_array_count(channel, preamble_array_number(channel)) = 0;
                        % 将补零后的峰值位置放入对应信道的preamble队列中，以便后续计算cfo和winoff
                        preamble_array_bin(channel, preamble_array_number(channel)) = pos;
                        % 将补零后的峰值强度也放入队列中，一遍后续分类SFD和downchirp
                        preamble_array_amp(channel, preamble_array_number(channel)) = peak;
                        % 将detect队列中对应元素删除，和上面做法相同
                        detect_array(channel, pos_count:end) = [detect_array(channel, pos_count+1:end), 0];
                        detect_array_count(channel, pos_count:end) = [detect_array_count(channel, pos_count+1:end), 0];
                        detect_array_number(channel) = detect_array_number(channel) - 1;
                    end
                end
            end
            % 对于那些没有检测到有效信号的信道，直接将其检测队列置零
            no_active_channel = find(~ismember([1,2,3,4], active_channel));
            % 遍历每一个此次没有检测到有效信号的信道
            for channel_num = 1:size(no_active_channel,2)
                no_channel = no_active_channel(channel_num);
                % 对其detect队列置零
                detect_array(no_channel, :) = zeros(1, 40);
                detect_array_count(no_channel, :) = zeros(1, 40);
                detect_array_number(no_channel) = 0;
            end
        end
        % 对preamble队列操作，找到sfd对应的bin
        % 每滑动一个窗口，给所以已经记录的preamble队列中所有数加1
        for channel_count = 1:size(signalOut,1)
            num = preamble_array_number(channel_count);
            if num >= 1
                preamble_array_count(channel_count, 1:num) = preamble_array_count(channel_count, 1:num) + 1;
            end
            % 对preamble队列中的每个峰进行遍历，如果count数目等于5开始计算downchirp_bin并计算cfo和winoff
            for pos_count = 1:num
                if preamble_array_count(channel_count, pos_count) == 5
                    % 对7个preamble后的第五个窗口进行补零FFT
                    [samples_fft_merge] = get_zeropadding_fft(signalOut(channel_count, 1:dine), loraSet, upchirp);
                    % 从preamble队列中拿到幅值信息，以便从sfd的fft从取得最相似的峰
                    amp_ref = preamble_array_amp(channel_count, pos_count);
                    % 获得去掉旁瓣影响的前N个峰
                    [pos_array_padding, peak_array_padding] = fft_filterN_peak(loraSet, samples_fft_merge, 40);
                    % 获得peak_array_padding与amp_ref的差，找到最接近的峰值
                    dif_array = abs(peak_array_padding - amp_ref);
                    % 找到最近的峰值
                    [~, min_pos] = min(dif_array);
                    % 获得downchirp_bin
                    downchirp_bin = pos_array_padding(min_pos);
                    % 从preamble队列中获得upchirp_bin
                    upchirp_bin = preamble_array_bin(channel_count, pos_count);
                    % 根据downchirp和upchirp bin获取cfo和winoff
                    [cfo, windows_offset] = get_cfo_winoff_bybin(loraSet, upchirp_bin, downchirp_bin);
                    % 将记录到的winoff和cfo存入到channel_array中，用于后续读取跳信道矩阵策略
                    channel_array_num(channel_count) = channel_array_num(channel_count) + 1;
                    pos = channel_array_num(channel_count);
                    channel_array_count(channel_count, pos) = 0;
                    channel_array_cfo(channel_count, pos) = cfo;
                    channel_array_winoff(channel_count, pos) = windows_offset;
                    channel_array_amp(channel_count, pos) = amp_ref;
                    % 在preamble队列中remove
                    preamble_array_bin(channel_count, pos_count:end) = [preamble_array_bin(channel_count, pos_count+1:end), 0];
                    preamble_array_amp(channel_count, pos_count:end) = [preamble_array_amp(channel_count, pos_count+1:end), 0];
                    preamble_array_count(channel_count, pos_count:end) = [preamble_array_count(channel_count, pos_count+1:end), 0];
                    preamble_array_number(channel_count) = preamble_array_number(channel_count) - 1;
                end
            end
        end
        % 对channel队列进行处理，获得信道矩阵策略的bin
        for channel_count = 1:size(signalOut,1)
            % 每滑动一个窗口，给所以已经记录的preamble队列中所有数加1
            num = channel_array_num(channel_count);
            if num >= 1
                channel_array_count(channel_count, 1:num) = channel_array_count(channel_count, 1:num) + 1;
            end
            % 对channel队列中的每个峰进行遍历，如果count数目等于3开始计算信道策略对应的bin值
            for pos_count = 1:num
                % 当count==5时，表示到达了获取信道策略的downchirp的窗口，开始解码
                if channel_array_count(channel_count, pos_count) == 5
                    % 获得记录在channel_array队列中的windows_offset
                    win_off = channel_array_winoff(channel_count, pos_count);
                    % 获得记录在channel_array队列中的cfo
                    cfo = channel_array_cfo(channel_count, pos_count);
                    % 对win_off取整
                    win_off = round(win_off);
                    % 对齐窗口
                    signal_tmp = signalOut(channel_count, dine*0.25+win_off+1: win_off+dine*1.25);
                    % 获得对应decfo的upchirp
                    [~, d_upchirp_cfo] = rebuild_idealchirp_cfo(loraSet, cfo, 0);
                    % 获取峰值幅度的参考值
                    amp_ref = channel_array_amp(channel_count, pos_count);
                    % 获得跳信道矩阵
                    [channel_jump_array] = get_DownchirpSync_multi(loraSet, signal_tmp, d_upchirp_cfo, amp_ref);
                    % 将信息存入解码队列（demodulate_array）中
                    demodulate_array_num = demodulate_array_num + 1;
                    demodulate_array_channel(demodulate_array_num, 1:200) = channel_jump_array(1:200);
                    demodulate_array_cfo(demodulate_array_num) = cfo;
                    demodulate_array_winoff(demodulate_array_num) = win_off;
                    % 将channel队列中的信息remove
                    channel_array_cfo(channel_count, pos_count:end) = [channel_array_cfo(channel_count, pos_count+1:end), 0];
                    channel_array_winoff(channel_count, pos_count:end) = [channel_array_winoff(channel_count, pos_count+1:end), 0];
                    channel_array_amp(channel_count, pos_count:end) = [channel_array_amp(channel_count, pos_count+1:end), 0];
                    channel_array_count(channel_count, pos_count:end) = [channel_array_count(channel_count, pos_count+1:end), 0];
                    channel_array_num(channel_count) = channel_array_num(channel_count) - 1;
                end
            end
        end
        % 对demodulate队列进行遍历，每个窗口都进行解码
        % 如果demodulate队列中存在内容，则每经过一个窗口，计数值加1
        if demodulate_array_num >= 1
            demodulate_array_count(1:demodulate_array_num) = demodulate_array_count(1:demodulate_array_num) + 1;
        end
        % 遍历demodulate队列
        for pos_count = 1:demodulate_array_num
            % 到达需要再次对齐的窗口
            if demodulate_array_count(pos_count) == 2
                % 从demodulate_array_winoff获取windows_off
                win_off = demodulate_array_winoff(pos_count);
                % 对齐窗口
                signal_tmp = signalOut(:, dine*0.25+win_off+1: win_off+dine*1.25);
                % 从demodulate_array_cfo获取cfo
                cfo = demodulate_array_cfo(pos_count);
                % 获得对应decfo的upchirp
                [d_downchirp_cfo, d_upchirp_cfo] = rebuild_idealchirp_cfo(loraSet, cfo, 0);
                % 从demodulate_array_channel获得该信号的跳信道矩阵
                channel_jump_array = demodulate_array_channel(pos_count, :);
                % 根据第一个跳信道subchirp（bin为0）来重新对齐
                [time_off] = align_windows_bysubchirp(loraSet, signal_tmp, channel_jump_array, 4, d_downchirp_cfo);
                % 将timeoff记录下来
                demodulate_array_timeoff(pos_count) = time_off;
                % 在此窗口下解调第一个chirp，看对齐是否正确
                % 重新对齐
                signal_tmp = signalOut(:, dine*0.25+win_off-time_off+1: win_off+dine*1.25-time_off);
                % 获得此次解调的4个subchirp的跳信道矩阵
                count = demodulate_array_count(pos_count) - 1;
                jump_channel = channel_jump_array((count-1)*4+1 : count*4);
                % 解调获得bin
                [subchirp_bin] = demodulate_bybin(loraSet, signal_tmp, jump_channel, 4, d_downchirp_cfo);
                % 将解得的bin记录下来
                demodulate_array_bin(pos_count, (count-1)*4+1 : count*4) = subchirp_bin;
            end
            if demodulate_array_count(pos_count) > 2 && demodulate_array_count(pos_count) <= 17
                % 从demodulate_array_winoff获取windows_off
                win_off = demodulate_array_winoff(pos_count);
                % 从demodulate_array_cfo获取cfo
                cfo = demodulate_array_cfo(pos_count);
                % 获得对应decfo的upchirp
                [d_downchirp_cfo, d_upchirp_cfo] = rebuild_idealchirp_cfo(loraSet, cfo, 0);
                % 从demodulate_array_channel获得该信号的跳信道矩阵
                channel_jump_array = demodulate_array_channel(pos_count, :);
                % 从demodulate_array_timeoff中获取time_offset
                time_off = demodulate_array_timeoff(pos_count);
                % 窗口对齐
                signal_tmp = signalOut(:, dine*0.25+win_off-time_off+1: win_off+dine*1.25-time_off);
                % 获得此次解调的4个subchirp的跳信道矩阵
                count = demodulate_array_count(pos_count) - 1;
                jump_channel = channel_jump_array((count-1)*4+1 : count*4);
                % 解调获得bin
                [subchirp_bin] = demodulate_bybin(loraSet, signal_tmp, jump_channel, 4, d_downchirp_cfo);
                % 将解得的bin记录下来
                demodulate_array_bin(pos_count, (count-1)*4+1 : count*4) = subchirp_bin;
            end
            % 当解调完16个bin，将其解调的bin记录下来，并从demodulate队列中删除
            if demodulate_array_count(pos_count) == 17
                % 记录bin
                demodulated_array_num = demodulated_array_num + 1;
                demodulated_array(demodulated_array_num, 1:16*4) = demodulate_array_bin(pos_count, 1:16*4);
                % 从demodulate队列中remove
                demodulate_array_bin(pos_count:end, :) = [demodulate_array_bin(pos_count+1:end, :); zeros(1, size(demodulate_array_bin, 2))];
                demodulate_array_channel(pos_count:end, :) = [demodulate_array_channel(pos_count+1:end, :); zeros(1, size(demodulate_array_channel, 2))];
                demodulate_array_cfo(pos_count : end) = [demodulate_array_cfo(pos_count+1 : end), 0];
                demodulate_array_winoff(pos_count : end) = [demodulate_array_winoff(pos_count+1 : end), 0];
                demodulate_array_count(pos_count : end) = [demodulate_array_count(pos_count+1 : end), 0];
                demodulate_array_timeoff(pos_count : end) = [demodulate_array_timeoff(pos_count+1 : end), 0];
                demodulate_array_num = demodulate_array_num - 1;
            end
        end
    end
end

toc;
fclose all;