% 画出所选文件的STFT，并给出bin值
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
leakage_width_array = [0.05,0.01,0.015,0.01,0.01,0.01];
lora_set.filter_num = 2;
lora_set.leakage_width1 = leakage_width_array(lora_set.sf-6);
lora_set.leakage_width2 = 1-lora_set.leakage_width1;
lora_set.channel = 4;
lora_set.channel_choice = [0,1,2,3];
lora_set.channel_choice_num = 4;
lora_set.Preamble_channel = 1;
Pkg_Verification = load(Verification_path)';
lora_set.pkg_vertification = Pkg_Verification;
[d_downchirp, d_upchirp] = build_idealchirp(lora_set, -187.5e3); % build idealchirp
% [d_downchirp_cfo, d_upchirp_cfo] = rebuild_idealchirp_cfo(lora_set, 1e5, 5e5);
% stft(d_downchirp_cfo(1:1*lora_set.dine), lora_set.sample_rate, 'Window',rectwin(256),'OverlapLength',128,'FFTLength',lora_set.fft_x);

Samples_Path = 'E:\share\';            % 设置采样值文件所在路径
File_1 = dir(fullfile(Samples_Path,'*.sigmf-data'));          % 读取文件夹下满足规则的采样值文件1

dine = lora_set.dine;
fft_x = lora_set.fft_x;
factor = dine/fft_x;
% for i = 1:length(File_1)
for i = 1:1

    if mod(i,100) == 0           % 每循环100次，输出当前进度
        fprintf("The time is %d\n",i);
    end
    % 读取数据
    File1_Path = strcat(Samples_Path, File_1(i).name);
    fid_1=fopen(File1_Path,'rb');
    [A_1]=fread(fid_1,'float32')';
    A_1_length = size(A_1,2);
    G0 = A_1(1:2:A_1_length-1) + A_1(2:2:A_1_length)*1i;
    % 找到preamble开始位置
    [Preamble_start_pos] = detect_preamble(G0, lora_set, d_downchirp); 
    detect_flag = 1;
    if Preamble_start_pos == 999
        detect_flag = 0;
    elseif Preamble_start_pos ~= 1
        G0 = circshift(G0, -(Preamble_start_pos-1) * lora_set.dine);
    end
    % 调整CFO和timeoff
    [cfo, windows_offset] = get_cfo_winoff(G0, lora_set, d_downchirp, d_upchirp, factor);
    G0 = circshift(G0, -round(windows_offset));
    figure(1);
    stft(G0(1:18.25*lora_set.dine), lora_set.sample_rate, 'Window',rectwin(256),'OverlapLength',128,'FFTLength',lora_set.fft_x);
    
    [downchirp_3175_cfo, upchirp_3175_cfo] = rebuild_idealchirp_cfo(lora_set, cfo, -187.5e3);
    [downchirp_3375_cfo, upchirp_3375_cfo] = rebuild_idealchirp_cfo(lora_set, cfo, -387.5e3);
    [downchirp_3575_cfo, upchirp_3575_cfo] = rebuild_idealchirp_cfo(lora_set, cfo, -587.5e3);
    [downchirp_3775_cfo, upchirp_3775_cfo] = rebuild_idealchirp_cfo(lora_set, cfo, -787.5e3);
%     G_tmp = [upchirp_3175_cfo,upchirp_3375_cfo,upchirp_3575_cfo,upchirp_3775_cfo];
%     G_tmp = upchirp_3775_cfo + upchirp_4175_cfo + upchirp_3175_cfo + circshift(upchirp_3575_cfo, length(upchirp_3575_cfo)/2);
%     G_tmp = G_tmp(1:16:end);
%     figure(2);
%     stft(G_tmp, lora_set.sample_rate, 'Window',rectwin(256/16),'OverlapLength',128/16,'FFTLength',lora_set.fft_x/16);
%     stft(G_tmp, lora_set.sample_rate, 'Window',rectwin(256),'OverlapLength',128,'FFTLength',lora_set.fft_x);
    %     figure(10);
%     stft(G0(12.25*dine+1:18.25*lora_set.dine), lora_set.sample_rate, 'Window',rectwin(256),'OverlapLength',128,'FFTLength',lora_set.fft_x);
% 
%     [downchirp_3175_cfo, upchirp_3175_cfo] = rebuild_idealchirp_cfo(lora_set, cfo, 5e5);
%     [downchirp_3375_cfo, upchirp_3375_cfo] = rebuild_idealchirp_cfo(lora_set, cfo, 3e5);
%     [downchirp_3575_cfo, upchirp_3575_cfo] = rebuild_idealchirp_cfo(lora_set, cfo, 1e5);
%     [downchirp_3775_cfo, upchirp_3775_cfo] = rebuild_idealchirp_cfo(lora_set, cfo, -1e5);
%     [downchirp_3975_cfo, upchirp_3975_cfo] = rebuild_idealchirp_cfo(lora_set, cfo, -3e5);
%     [downchirp_4175_cfo, upchirp_4175_cfo] = rebuild_idealchirp_cfo(lora_set, cfo, -5e5);

    [singal_out] = divide_channel(lora_set, G0);
    figure(2);
    stft(singal_out(1, 1:18.25*lora_set.dine), lora_set.sample_rate, 'Window',rectwin(256),'OverlapLength',128,'FFTLength',lora_set.fft_x);
    figure(3);
    stft(singal_out(2, 1:18.25*lora_set.dine), lora_set.sample_rate, 'Window',rectwin(256),'OverlapLength',128,'FFTLength',lora_set.fft_x);
    figure(4);
    stft(singal_out(3, 1:18.25*lora_set.dine), lora_set.sample_rate, 'Window',rectwin(256),'OverlapLength',128,'FFTLength',lora_set.fft_x);
    figure(5);
    stft(singal_out(4, 1:18.25*lora_set.dine), lora_set.sample_rate, 'Window',rectwin(256),'OverlapLength',128,'FFTLength',lora_set.fft_x);

end

toc;
fclose all;