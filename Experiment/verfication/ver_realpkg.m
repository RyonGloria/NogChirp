% 将所选文件夹中的所有*.sigmf-data文件对齐后写入到目的文件夹
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

Config_Path = '.\config\';                                      % 设置配置文件所在路径
% Verification_path = strcat(Config_Path,'helloworld.txt');       % bin值验证文件
Setting_File = dir(fullfile(Config_Path,'sf10_BW125.json'));    % 配置文件
Setting_File_Path = strcat(Config_Path, Setting_File.name);
Setting_file = fopen(Setting_File_Path,'r');
setting = jsondecode(fscanf(Setting_file,'%s'));                % 解析json格式变量
lora_set.bw = setting.captures.lora_bw;                         % 设置接收数据包的lora_set.bw
lora_set.sf = setting.captures.lora_sf;                         % 设置接收数据包的lora_set.sf
lora_set.sample_rate = setting.global.core_sample_rate;         % 设置接收数据包的采样率
% lora_set.sample_rate = 500000;
lora_set.Pkg_length = setting.captures.lora_pkg_length;         % 设置接收数据包的长度
lora_set.dine = lora_set.sample_rate*bitshift(1,lora_set.sf)/lora_set.bw;    % 根据lora_set.sf和lora_set.bw计算出一个chirp包含的采样点个数
lora_set.fft_x = 2^lora_set.sf;                               % 根据lora_set.dine计算出包含lora_set.bw所需的FFT点数
lora_set.Preamble_length = 8;
leakage_width_array = [0.05,0.01,0.015,0.01];
lora_set.filter_num = 2;
lora_set.leakage_width1 = leakage_width_array(lora_set.sf-6);
lora_set.leakage_width2 = 1-lora_set.leakage_width1;
lora_set.channel = 1;
lora_set.channel_choice = [0];
lora_set.channel_choice_num = 1;
lora_set.Preamble_channel = 0;

bin_string = '257';
preamble_string = string(lora_set.Preamble_channel);
channel_string = string(lora_set.channel_choice_num);

% Samples_Path = strcat('D:\CHchirp_samples\align_windows\stable_binchannel\channel',channel_string,'\preamble',preamble_string,'\bin',bin_string,'\');
Samples_Path = strcat('D:\CHchirp_samples\align_windows\rand_binchannel\channel',channel_string,'\');
% Samples_Path = strcat('D:\CHchirp_samples\rand_binchannel\channel',channel_string,'\');
% Verification_path = strcat('D:\CHchirp_samples\stable_binchannel_GT\channel',channel_string,'\preamble',preamble_string,'\bin',bin_string,'\');
Verification_path = strcat('D:\CHchirp_samples\rand_binchannel_GT\channel',channel_string,'\');
File_1 = dir(fullfile(Samples_Path,'*.sigmf-data'));          % 读取文件夹下满足规则的采样值文件1
VerFile = dir(fullfile(Verification_path,'bin_*.txt'));  
VerchannelFile = dir(fullfile(Verification_path,'channel_*.txt'));
datenum_array = zeros(1, length(File_1));
ver_date_array = zeros(1, length(VerFile));
vercha_date_array = zeros(1, length(VerchannelFile));
for i = 1:length(File_1)
    datenum_array(i) = File_1(i).datenum;
    ver_date_array(i) = VerFile(i).datenum;
    vercha_date_array(i) = VerchannelFile(i).datenum;
end
[~, date_sort_index] = sort(datenum_array);
[~, ver_date_sort_index] = sort(ver_date_array);
[~, vercha_date_sort_index] = sort(vercha_date_array);
error_count = 0;
for i = 1:length(File_1)
%     if mod(i,100) == 0           % 每循环100次，输出当前进度
%         fprintf("The time is %d\n",i);
%     end
    % 读取GT
    Verbin_file_path = strcat(Verification_path, VerFile(ver_date_sort_index(i)).name);
    Pkg_Verification = load(Verbin_file_path)';
    Verchannel_file_path = strcat(Verification_path, VerchannelFile(vercha_date_sort_index(i)).name);
    Channel_Verification = load(Verchannel_file_path);
    % 读取文件
    File1_Path = strcat(Samples_Path, File_1(date_sort_index(i)).name);
    fid_1=fopen(File1_Path,'rb');
    [A_1]=fread(fid_1,'float32')';
    fclose all;
    A_1_length = size(A_1,2);
    G0 = A_1(1:2:A_1_length-1) + A_1(2:2:A_1_length)*1i;
    % 抽取成采样率为500KHz的信号
%     G0 = G0(1:4:end-3);
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
    % 分别解调preamble和payload，获取bin值
    [max_pos_array] = demodulate_pre(G0, lora_set, d_upchirp_cfo, d_downchirp_cfo);
    [bin_array, channel_choice_array] = get_channel_bin(G0, lora_set, cfo, 20);
    % [max_pos_array] = demodulate_pre(G0, lora_set, d_upchirp_cfo, d_downchirp_cfo);
    % Pkg_Verification = [1,1,1,1,1,1,1,1,9,17];
    bin_true_num = sum(bin_array == Pkg_Verification);
    channel_true = sum(channel_choice_array == Channel_Verification(:,1),"all");
    if channel_true == 20 && bin_true_num == 20
        
    else
        disp(File1_Path);

        delete(File1_Path);
        delete(Verbin_file_path);
        delete(Verchannel_file_path);
    end
    
    fclose all;
end
% disp(Samples_Path);
% disp(write_dir);
disp(error_count);

toc;
fclose all;