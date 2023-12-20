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
% lora_set.sample_rate = setting.global.core_sample_rate;         % 设置接收数据包的采样率
lora_set.sample_rate = 500000;
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

bin_string = '385';
preamble_string = string(lora_set.Preamble_channel);
channel_string = string(lora_set.channel_choice_num);

Samples_Path = strcat('D:\CHchirp_samples\align_windows\rand_binchannel\channel',channel_string,'\preamble',preamble_string,'\');
Verification_path = strcat('D:\CHchirp_samples\rand_binchannel_GT\channel',channel_string,'\preamble',preamble_string,'\bin',bin_string,'\');
File_1 = dir(fullfile(Samples_Path,'*.sigmf-data'));          % 读取文件夹下满足规则的采样值文件1
VerFile = dir(fullfile(Verification_path,'bin_*.txt'));  
VerchannelFile = dir(fullfile(Verification_path,'channel_*.txt'));
datenum_array = zeros(1, length(File_1));
ver_date_array = zeros(1, length(VerFile));
vercha_date_array = zeros(1, length(VerchannelFile));
for i = 1:length(File_1)
    datenum_array(i) = File_1(i).datenum;
%     ver_date_array(i) = VerFile(i).datenum;
%     vercha_date_array(i) = VerchannelFile(i).datenum;
end
[~, date_sort_index] = sort(datenum_array);
% [~, ver_date_sort_index] = sort(ver_date_array);
% [~, vercha_date_sort_index] = sort(vercha_date_array);
error_count = 0;
for i = 1:1
    % 读取文件
    File1_Path = strcat(Samples_Path, File_1(date_sort_index(i)).name);
    fid_1=fopen(File1_Path,'rb');
    [A_1]=fread(fid_1,'float32')';
    A_1_length = size(A_1,2);
    G0 = A_1(1:2:A_1_length-1) + A_1(2:2:A_1_length)*1i;
    % 抽取成采样率为500KHz的信号
    G0 = G0(1:4:end-3);
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
    [cfo, windows_offset] = get_cfo_winoff(G0, lora_set, d_downchirp, d_upchirp, 4);
    G0 = circshift(G0, -round(windows_offset));
    [d_downchirp_cfo, d_upchirp_cfo] = rebuild_idealchirp_cfo_pre(lora_set, cfo);
    [d_downchirp_cfo_pay] = rebuild_idealchirp_cfo_pay(lora_set, cfo);
    % 分别解调preamble和payload，获取bin值
    payload_length = 1;
    dine = lora_set.dine;
    fft_x = lora_set.fft_x;
    Ch_num = lora_set.channel_choice_num;
    start_pos = lora_set.Preamble_length+2+2.25;
    G0_pay = G0(start_pos*dine+1:(start_pos+payload_length)*dine);
    % 生成所有信道的downchirp
    samples = G0_pay(1:dine);
    lora_set.channel_choice = zeros(1, Ch_num);
    downchirp_tmp = rebuild_idealchirp_cfo_pay(lora_set, cfo);
    dechirp = samples.*downchirp_tmp;
    [~,bin_tmp] = max(abs(fft(dechirp)));
    bin = mod(bin_tmp, fft_x);
    
    downchirp_array = zeros(Ch_num, dine);
    for i = 1:Ch_num
%         lora_set.channel_choice = repmat(i-1, 1, Ch_num);
        lora_set.channel_choice = [2,3,0,0];
        downchirp_tmp = rebuild_idealchirp_cfo_pay(lora_set, cfo);
        [max_pos_array] = demodulate_pay(G0, lora_set, downchirp_tmp, 1);
%         dechirp = samples.*downchirp_array(i, :);
%         subplot(4,1,i);
%         stft(dechirp, lora_set.sample_rate, 'Window',rectwin(128),'OverlapLength',0,'FFTLength',2048);
    end
    % 对每个chirp遍历所有信道可能，找到正确的bin和信道策略
    fclose all;
end


toc;
fclose all;