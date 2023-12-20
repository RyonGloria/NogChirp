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
lora_set.channel = 2;
lora_set.channel_choice = [0,1];
lora_set.channel_choice_num = 2;
lora_set.Preamble_channel = 0;

% 生成理想仿真信号
payload_bin = [0,1,2,3,4,5,6,7,8,9,10,1024,1023,1022,1021,1020];
[d_downchirp, d_upchirp] = build_idealchirp_pre(lora_set);
[G0] = build_sim_lorapkg(lora_set, payload_bin);
% 检测preamble
[Preamble_start_pos] = detect_preamble(G0, lora_set, d_downchirp);
if Preamble_start_pos ~= 1
    G0 = circshift(G0, -(Preamble_start_pos-1) * lora_set.dine);
end

stft(G0(1:20*lora_set.dine), lora_set.sample_rate, 'Window',rectwin(128),'OverlapLength',0,'FFTLength',2048)
% 对齐窗口和cfo
% [cfo, windows_offset] = get_cfo_winoff(G0, lora_set, d_downchirp, d_upchirp, 8);
% G0 = circshift(G0, -round(windows_offset));
% [d_downchirp_cfo, d_upchirp_cfo] = rebuild_idealchirp_cfo_pre(lora_set, cfo);
% [d_downchirp_cfo_pay] = rebuild_idealchirp_cfo_pay(lora_set, cfo);
% % 分别解调preamble，获得其bin值
% [max_pos_array] = demodulate_pre(G0, lora_set, d_upchirp_cfo, d_downchirp_cfo);
% disp(max_pos_array');
% [max_pos_array] = demodulate_pay(G0, lora_set, d_downchirp_cfo_pay, 16);
% disp(max_pos_array');