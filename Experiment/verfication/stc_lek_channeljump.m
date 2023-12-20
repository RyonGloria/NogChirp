% 统计4信道情况下，跳信道信号的平均能量统计
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

Config_Path = '.\config\';                                      % 设置配置文件所在路径
Verification_path = strcat(Config_Path,'helloworld.txt');       % bin值验证文件
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
lora_set.channel = 4;
lora_set.channel_choice = [0,1,2,3];
lora_set.channel_choice_num = 4;
lora_set.Preamble_channel = 0;
Pkg_Verification = load(Verification_path)';
lora_set.pkg_vertification = Pkg_Verification;

Samples_Path = 'D:\CHchirp_samples\align_windows\USRP\4channel0\';                    % 设置采样值文件所在路径
File_1 = dir(fullfile(Samples_Path,'4channel0000*.sigmf-data'));          % 读取文件夹下满足规则的采样值文件1

peak_energy_all = zeros(1,16);
peak_energy_count = zeros(1,16);
channel_str = ['4channel0000*.sigmf-data';'4channel0011*.sigmf-data';'4channel0022*.sigmf-data';'4channel0033*.sigmf-data'];
channel_str = [channel_str; '4channel1100*.sigmf-data';'4channel1111*.sigmf-data';'4channel1122*.sigmf-data';'4channel1133*.sigmf-data'];
channel_str = [channel_str; '4channel2200*.sigmf-data';'4channel2211*.sigmf-data';'4channel2222*.sigmf-data';'4channel2233*.sigmf-data'];
channel_str = [channel_str; '4channel3300*.sigmf-data';'4channel3311*.sigmf-data';'4channel3322*.sigmf-data';'4channel3333*.sigmf-data'];
for channel_count = 1:16
    Samples_Path = strcat('D:\CHchirp_samples\align_windows\USRP\4channel', string(floor((channel_count-1)/4)), '\');
    File_1 = dir(fullfile(Samples_Path, channel_str(channel_count,:)));
    for i = 1:length(File_1)
        % 根据文件设置信道选择方案和preamble信道
        tmp = File_1(i).name;
        lora_set.channel_choice = [str2double(tmp(9)), str2double(tmp(10)), str2double(tmp(11)), str2double(tmp(12))];
        lora_set.Preamble_channel = str2double(tmp(9));
    
        % 读取文件
        File1_Path = strcat(Samples_Path, File_1(i).name);
        fid_1=fopen(File1_Path,'rb');
        [A_1]=fread(fid_1,'float32')';
        A_1_length = size(A_1,2);
        G0 = A_1(1:2:A_1_length-1) + A_1(2:2:A_1_length)*1i;
        % build ideachirp
        [d_downchirp, d_upchirp] = build_idealchirp_pre(lora_set);   
        % 检测preamble
        [Preamble_start_pos] = detect_preamble(G0, lora_set, d_downchirp); 
        detect_flag = 1;
        if Preamble_start_pos == 999
            detect_flag = 0;
        elseif Preamble_start_pos ~= 1
            G0 = circshift(G0, -(Preamble_start_pos-1) * lora_set.dine);
        end
        % 计算cfo和窗口对齐
        [cfo, windows_offset] = get_cfo_winoff(G0, lora_set, d_downchirp, d_upchirp, 16);
        G0 = circshift(G0, -round(windows_offset));
        [d_downchirp_cfo, d_upchirp_cfo] = rebuild_idealchirp_cfo_pre(lora_set, cfo);
        [d_downchirp_cfo_pay] = rebuild_idealchirp_cfo_pay(lora_set, cfo);
    
        % 解调payload，获取bin值和峰值能量
        dine = lora_set.dine;
        fft_x = lora_set.fft_x;
        payload_length = 23;
        start_pos = lora_set.Preamble_length+2+2.25;
        samples = reshape(G0(start_pos*dine+1:(start_pos+payload_length)*dine), [dine, payload_length]).';
    
        samples_dechirp = samples .* d_downchirp_cfo_pay;
        samples_dechirp(:, 1:30) = 0;
        samples_dechirp(:, dine-30:dine) = 0;
        samples_fft = abs(fft(samples_dechirp,dine,2));
        samples_fft_merge = samples_fft(:,1:fft_x) + samples_fft(:,dine-fft_x+1:dine);
        [peak_energy, max_pos_array] = max(samples_fft_merge, [], 2);
        true_num = sum((max_pos_array-1)' == lora_set.pkg_vertification);
    
        if detect_flag == 1 && true_num == 23 % 作为统计正确的标志
            peak_energy_count(channel_count) = peak_energy_count(channel_count) + 23;
            peak_energy_all(channel_count) = peak_energy_all(channel_count) + sum(peak_energy);
        else
            disp(File1_Path);
        end
        
        fclose all;
    end
end
% 计算结果
peak_energy_real = peak_energy_all ./ peak_energy_count;
peak_energy_real = reshape(peak_energy_real, [4,4]).';
% 生成对比能量
peak_energy_sim = zeros(4,4);
for col = 1:4
    for row = 1:4
        peak_energy_sim(col, row) = 0.5*peak_energy_real(col,col) + 0.5*peak_energy_real(row,row);
    end
end

subplot(4,1,1);
plot(1:4, peak_energy_real(1,:)./peak_energy_sim(1,:), '');
xticks([1 2 3 4]);
xticklabels({'0->0','0->1','0->2','0->3'});
ylabel('peak_energy');
subplot(4,1,2);
plot(1:4,peak_energy_real(2,:)./peak_energy_sim(2,:),'');
xticks([1 2 3 4]);
xticklabels({'1->0','1->1','1->2','1->3'});
ylabel('peak_energy');
subplot(4,1,3);
plot(1:4,peak_energy_real(3,:)./peak_energy_sim(3,:),'');
xticks([1 2 3 4]);
xticklabels({'2->0','2->1','2->2','2->3'});
ylabel('peak_energy');
subplot(4,1,4);
plot(1:4,peak_energy_real(4,:)./peak_energy_sim(4,:),'');
xticks([1 2 3 4]);
xticklabels({'3->0','3->1','3->2','3->3'});
ylabel('peak_energy');

toc;
fclose all;