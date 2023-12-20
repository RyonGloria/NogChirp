%% 统计所选文件夹下所有信号文件的FFT能量频率泄露
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

Config_Path = '.\config\';                                      % 设置配置文件所在路径
% Verification_path = strcat(Config_Path,'sf10_BW500.txt');       % bin值验证文件
Setting_File = dir(fullfile(Config_Path,'sf10_BW125.json'));    % 配置文件
Setting_File_Path = strcat(Config_Path, Setting_File.name);
Setting_file = fopen(Setting_File_Path,'r');
setting = jsondecode(fscanf(Setting_file,'%s'));                % 解析json格式变量
lora_set.bw = setting.captures.lora_bw;                         % 设置接收数据包的lora_set.bw
lora_set.sf = setting.captures.lora_sf;                         % 设置接收数据包的lora_set.sf
lora_set.sample_rate = setting.global.core_sample_rate;         % 设置接收数据包的采样率
lora_set.sample_rate = 2e6;

lora_set.Pkg_length = setting.captures.lora_pkg_length;         % 设置接收数据包的长度
lora_set.dine = lora_set.sample_rate*bitshift(1,lora_set.sf)/lora_set.bw;    % 根据lora_set.sf和lora_set.bw计算出一个chirp包含的采样点个数
lora_set.fft_x = 2^lora_set.sf;                               % 根据lora_set.dine计算出包含lora_set.bw所需的FFT点数
lora_set.Preamble_length = 8;
leakage_width_array = [0.05,0.01,0.015,0.01];
lora_set.filter_num = 2;
lora_set.leakage_width1 = leakage_width_array(lora_set.sf-6);
lora_set.leakage_width2 = 1-lora_set.leakage_width1;

Samples_Path = 'E:\';   % 设置采样值文件所在路径
File_1 = dir(fullfile(Samples_Path,'*.sigmf-data'));             % 读取文件夹下满足规则的采样值文件
% [d_downchirp, d_upchirp] = build_idealchirp(lora_set); % build idealchirp
[d_downchirp, d_upchirp] = build_idealchirp(lora_set, 5e5); % build idealchirp
result_array = zeros(3, length(File_1)*100);    % 设置变量记录结果数据
bin_diff = zeros(0);
bin_peak = zeros(0);

dine = lora_set.dine;
fft_x = lora_set.fft_x;
for i = 1:length(File_1)
    if mod(i,100) == 0           % 每循环100次，输出当前进度
        fprintf("The time is %d\n",i);
    end
    % 读取文件
    File1_Path = strcat(Samples_Path, File_1(i).name);
    fid_1=fopen(File1_Path,'rb');
    [A_1]=fread(fid_1,'float32')';
    A_1_length = size(A_1,2);
    G0 = A_1(1:2:A_1_length-1) + A_1(2:2:A_1_length)*1i;
    %找到preamble开始位置
    [Preamble_start_pos] = detect_preamble(G0, lora_set, d_downchirp); 
    detect_flag = 1;
    if Preamble_start_pos == 999
        detect_flag = 0;
    elseif Preamble_start_pos ~= 1
        G0 = circshift(G0, -(Preamble_start_pos-1) * lora_set.dine);
    end
    % 计算cfo和time_offset
    [cfo, windows_offset] = get_cfo_winoff(G0, lora_set, d_downchirp, d_upchirp, 16);
    G0 = circshift(G0, -round(windows_offset));
%     [d_downchirp_cfo, d_upchirp_cfo] = rebuild_idealchirp_cfo(lora_set, cfo);
    [downchirp_3175_cfo, upchirp_3175_cfo] = rebuild_idealchirp_cfo(lora_set, cfo, 5e5);
    [downchirp_3375_cfo, upchirp_3375_cfo] = rebuild_idealchirp_cfo(lora_set, cfo, 3e5);
    [downchirp_3575_cfo, upchirp_3575_cfo] = rebuild_idealchirp_cfo(lora_set, cfo, 1e5);
    [downchirp_3775_cfo, upchirp_3775_cfo] = rebuild_idealchirp_cfo(lora_set, cfo, -1e5);
    [downchirp_3975_cfo, upchirp_3975_cfo] = rebuild_idealchirp_cfo(lora_set, cfo, -3e5);
    [downchirp_4175_cfo, upchirp_4175_cfo] = rebuild_idealchirp_cfo(lora_set, cfo, -5e5);
    downchirp_cfo = [downchirp_3175_cfo; downchirp_3375_cfo; downchirp_3575_cfo; downchirp_3775_cfo; downchirp_3975_cfo; downchirp_4175_cfo];
    upchirp_cfo = [upchirp_3175_cfo; upchirp_3375_cfo; upchirp_3575_cfo; upchirp_3775_cfo; upchirp_3975_cfo; upchirp_4175_cfo];
    % 根据相位精对齐
%     result = zeros(0);
%     result_bin = zeros(0);
%     G0_tmp = G0(20.75*dine+1 : 21.75*dine);
%     dechirp = G0_tmp .* d_downchirp_cfo;
%     dechirp_fft = fft(dechirp, lora_set.dine);
%     for ang = 0 : 2*pi/360 : 2*pi-2*pi/360
%         dechirp_fft_merge = dechirp_fft(1:fft_x)*exp(ang*1j) + dechirp_fft(fft_x+1:dine);
%         [tmp, tmp_pos] = max(dechirp_fft_merge);
%         result = [result, tmp];
%         result_bin = [result_bin, tmp_pos];
%     end
%     [~, tmp_pos] = max(result);
%     if tmp_pos < 60 || tmp_pos > 120
%         G0 = circshift(G0, 1);
%     end
    
    % 统计每个chirp的FFT能量频率泄露
    hop_count = 1;
    for chirp_count = 13:2:211
        G0_tmp_1 = G0((chirp_count-1)*dine+0.25*dine+1 : chirp_count*dine+0.25*dine);
        G0_tmp_2 = G0((chirp_count)*dine+0.25*dine+1 : (chirp_count+1)*dine+0.25*dine);
        % 与对应信道的idealchirp解调
        if mod(hop_count,2) == 1
            dechirp = [G0_tmp_1.*downchirp_cfo(1, :), G0_tmp_2.*downchirp_cfo(ceil(hop_count/2)+1, :)];   % 分别对两个chirp进行dechirp
        else
            dechirp = [G0_tmp_1.*downchirp_cfo(ceil(hop_count/2)+1, :), G0_tmp_2.*downchirp_cfo(1, :)];
        end
        
        % 划分四个窗口
        dechirp_1_1 = [zeros(1, 20), dechirp(21:dine)];
        dechirp_1_2 = [dechirp(dine+1:2*dine-20), zeros(1, 20)];
        dechirp_2 = [zeros(1, 10), dechirp(11:dine-10), zeros(1, 10)];
        dechirp_3 = [zeros(1, 10), dechirp(dine+11:2*dine-10), zeros(1, 10)];
        % 分别对这四个窗口进行fft
        dechirp_1_1_fft = abs(fft(dechirp_1_1));
        dechirp_1_2_fft = abs(fft(dechirp_1_2));
        dechirp_2_fft = abs(fft(dechirp_2));
        dechirp_3_fft = abs(fft(dechirp_3));
        % fft_merge
        dechirp_1_1_fft = dechirp_1_1_fft(1:fft_x) + dechirp_1_1_fft(dine-fft_x+1:dine);
        dechirp_1_2_fft = dechirp_1_2_fft(1:fft_x) + dechirp_1_2_fft(dine-fft_x+1:dine);
        dechirp_2_fft = dechirp_2_fft(1:fft_x) + dechirp_2_fft(dine-fft_x+1:dine);
        dechirp_3_fft = dechirp_3_fft(1:fft_x) + dechirp_3_fft(dine-fft_x+1:dine);
        % 获得两个chirp的bin值
        [max_peak1, max_bin_1] = max(dechirp_2_fft);
        [max_peak2, max_bin_2] = max(dechirp_3_fft);
        % 找到他们对应中间窗口的峰
        condition_1 = abs(max_bin_1-[1:fft_x]) < fft_x*0.01;   % 找到fft_x*leakage_width1范围内的旁瓣
        condition_2 = abs(max_bin_1-[1:fft_x]) > fft_x*0.99;
        sidelobe_index = (condition_1 | condition_2);
        dechirp_1_max_left = max(dechirp_1_1_fft(sidelobe_index));
        condition_1 = abs(max_bin_2-[1:fft_x]) < fft_x*0.01;   % 找到fft_x*leakage_width1范围内的旁瓣
        condition_2 = abs(max_bin_2-[1:fft_x]) > fft_x*0.99;
        sidelobe_index = (condition_1 | condition_2);
        dechirp_1_max_right = max(dechirp_1_2_fft(sidelobe_index));
        % 计算峰值能量差
        diff_1 = (max_peak1 - dechirp_1_max_left)/max_peak1;
        diff_2 = (max_peak2 - dechirp_1_max_right)/max_peak2;
        % 计算两个的bin值差
        if mod(hop_count,2) == 1
            if max_bin_1 == 1
                fre_dif = 5e4 + (1.75e5)*(ceil(hop_count/2)-1) + lora_set.bw/fft_x * max_bin_2;
            else
                fre_dif = 5e4 + (1.75e5)*(ceil(hop_count/2)-1) + lora_set.bw/fft_x * (fft_x - max_bin_1) + lora_set.bw/fft_x * max_bin_2;
            end
            result_array(:, (i-1)*100 + ceil((chirp_count-12)/2) ) = [fre_dif; diff_1; diff_2];
        else
            if max_bin_1 == 1
                fre_dif = -5e4 - (1.75e5)*(ceil(hop_count/2)-1) - lora_set.bw - lora_set.bw/fft_x * (fft_x - max_bin_2);
            else
                fre_dif = -5e4 - (1.75e5)*(ceil(hop_count/2)-1) - lora_set.bw/fft_x * max_bin_1 - lora_set.bw/fft_x * (fft_x - max_bin_2);
            end
            result_array(:, (i-1)*100 + ceil((chirp_count-12)/2) ) = [fre_dif; diff_1; diff_2];
        end
        
        % 记录hop跳的次数
        hop_count = hop_count + 1;
        if(hop_count > 10)
            hop_count = 1;
        end
    end
%     for chirp_count = 16:130
%         G0_tmp_1 = G0((chirp_count-1)*dine+0.25*dine+1 : chirp_count*dine+0.25*dine);
%         G0_tmp_2 = G0((chirp_count)*dine+0.25*dine+1 : (chirp_count+1)*dine+0.25*dine);
%         dechirp = [G0_tmp_1.*d_downchirp_cfo, G0_tmp_2.*d_downchirp_cfo];   % 分别对两个chirp进行dechirp
%         
%         % 划分四个窗口
%         dechirp_1_1 = [dechirp(dine-2027:dine), zeros(1, 20)];
%         dechirp_1_2 = [dechirp(dine+1:dine+2028), zeros(1, 20)];
%         dechirp_2 = [dechirp(dine-2037:dine-10), zeros(1, 20)];
%         dechirp_3 = [dechirp(dine+11:dine+2038), zeros(1, 20)];
%         % 分别对这四个窗口进行fft
%         dechirp_1_1_fft = abs(fft(dechirp_1_1));
%         dechirp_1_2_fft = abs(fft(dechirp_1_2));
%         dechirp_2_fft = abs(fft(dechirp_2));
%         dechirp_3_fft = abs(fft(dechirp_3));
%         % fft_merge
%         dechirp_1_1_fft = dechirp_1_1_fft(1:fft_x) + dechirp_1_1_fft(dine-fft_x+1:dine);
%         dechirp_1_2_fft = dechirp_1_2_fft(1:fft_x) + dechirp_1_2_fft(dine-fft_x+1:dine);
%         dechirp_2_fft = dechirp_2_fft(1:fft_x) + dechirp_2_fft(dine-fft_x+1:dine);
%         dechirp_3_fft = dechirp_3_fft(1:fft_x) + dechirp_3_fft(dine-fft_x+1:dine);
%         % 获得两个chirp的bin值
%         [max_peak1, max_bin_1] = max(dechirp_2_fft);
%         [max_peak2, max_bin_2] = max(dechirp_3_fft);
%         % 找到他们对应中间窗口的峰
%         condition_1 = abs(max_bin_1-[1:fft_x]) < fft_x*0.01;   % 找到fft_x*leakage_width1范围内的旁瓣
%         condition_2 = abs(max_bin_1-[1:fft_x]) > fft_x*0.99;
%         sidelobe_index = (condition_1 | condition_2);
%         dechirp_1_max_left = max(dechirp_1_1_fft(sidelobe_index));
%         condition_1 = abs(max_bin_2-[1:fft_x]) < fft_x*0.01;   % 找到fft_x*leakage_width1范围内的旁瓣
%         condition_2 = abs(max_bin_2-[1:fft_x]) > fft_x*0.99;
%         sidelobe_index = (condition_1 | condition_2);
%         dechirp_1_max_right = max(dechirp_1_2_fft(sidelobe_index));
%         % 计算峰值能量差
%         diff_1 = (max_peak1 - dechirp_1_max_left)/max_peak1;
%         diff_2 = (max_peak2 - dechirp_1_max_right)/max_peak2;
%         % 计算两个的bin值差
%         fre_dif = (max_bin_1 - max_bin_2);
%         result_array(:, (i-1)*(130-16+1)+(chirp_count-15)) = [fre_dif; diff_1; diff_2];
%     end
    fclose all;
end

toc;
fclose all;