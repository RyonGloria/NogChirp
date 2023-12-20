% 画出所选文件的fft图
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

Config_Path = '.\config\';                                       % 设置配置文件所在路径
% Verification_path = strcat(Config_Path,'node1_sf10.txt');        % bin值验证文件
Setting_File = dir(fullfile(Config_Path,'sf7_BW125.json'));     % 配置文件
Setting_File_Path = strcat(Config_Path, Setting_File.name);
Setting_file = fopen(Setting_File_Path,'r');
setting = jsondecode(fscanf(Setting_file,'%s'));                % 解析json格式变量
lora_set.bw = setting.captures.lora_bw;                         % 设置接收数据包的lora_set.bw
lora_set.sf = setting.captures.lora_sf;                                  % 设置接收数据包的lora_set.sf
lora_set.sample_rate = setting.global.core_sample_rate;                  % 设置接收数据包的采样率
lora_set.Pkg_length = setting.captures.lora_pkg_length;                  % 设置接收数据包的长度
lora_set.dine = 1000000*bitshift(1,lora_set.sf)/lora_set.bw;                               % 根据lora_set.sf和lora_set.bw计算出一个chirp包含的采样点个数
lora_set.fft_x = lora_set.dine/2;                                                 % 根据lora_set.dine计算出包含lora_set.bw所需的FFT点数
lora_set.Preamble_length = 8;
leakage_width_array = [0.05,0.02,0.02,0.02];
lora_set.filter_num = 2;
lora_set.leakage_width1 = leakage_width_array(lora_set.sf-6);
lora_set.leakage_width2 = 1-lora_set.leakage_width1;

% 读取特定文件
fid_1=fopen('D:\CHchirp_samples\USRP\T21^%53^%15_SF10_BW500000_R18727.sigmf-data','rb');
[A_1]=fread(fid_1,'float32')';
A_1_length = size(A_1,2);
G0 = A_1(1:2:A_1_length-1) + A_1(2:2:A_1_length)*1i;            % 将float数组转换成复数数组
[d_downchirp, d_upchirp] = build_idealchirp(lora_set);          % 生成idealchirp
% 对齐并计算cfo
[cfo, windows_offset] = get_cfo_winoff(G0, lora_set, d_downchirp, d_upchirp, 2);
G0 = circshift(G0, -round(windows_offset));
[d_downchirp_cfo, d_upchirp_cfo] = rebuild_idealchirp_cfo(lora_set, cfo);

% 获取其中一个窗口的信号
dine = lora_set.dine;
fft_x = lora_set.fft_x;
G0_tmp = G0(20.75*dine+1 : 21.75*dine);
% time_plot(G0_tmp, lora_set, d_downchirp_cfo, 30);
% tmp = circshift(d_upchirp, 0.25*dine);
% tmp2 = circshift(d_upchirp, 0.75*dine);
% G0_tmp = [tmp(dine/2+1:end), tmp2(1:dine/2)];
% dechirp = G0_tmp .* d_downchirp;
dechirp = G0_tmp .* d_downchirp_cfo;
dechirp_fft = fft(dechirp, lora_set.dine);

% 通过相位精对齐
result = zeros(0);
result_bin = zeros(0);
for ang = 0 : 2*pi/360 : 2*pi-2*pi/360
    dechirp_fft_merge = dechirp_fft(1:fft_x)*exp(ang*1j) + dechirp_fft(fft_x+1:dine);
    [tmp, tmp_pos] = max(dechirp_fft_merge);
    result = [result, tmp];
    result_bin = [result_bin, tmp_pos];
end
[~, tmp_pos] = max(result);
if tmp_pos < 60 || tmp_pos > 120
    G0 = circshift(G0, 1);
end

% 裁剪三个窗口并画出他们的fft
G0_tmp_1 = G0(20.25*dine+1 : 21.25*dine);
G0_tmp_2 = G0(21.25*dine+1 : 22.25*dine);
dechirp = [G0_tmp_1.*d_downchirp_cfo, G0_tmp_2.*d_downchirp_cfo];
% dechirp_1 = dechirp(dine/2-128:dine/2+128);
% dechirp_2 = [dechirp(dine/2-256:dine/2-128), zeros(1,128)];
% dechirp_3 = [dechirp(dine/2+128:dine/2+256), zeros(1,128)];
dechirp_1_1 = [dechirp(dine-2027:dine), zeros(1, 20)];
dechirp_1_2 = [dechirp(dine+1:dine+2028), zeros(1, 20)];
dechirp_2 = [dechirp(dine-2037:dine-10), zeros(1, 20)];
dechirp_3 = [dechirp(dine+11:dine+2038), zeros(1, 20)];

dechirp_1_1_fft = abs(fft(dechirp_1_1));
dechirp_1_2_fft = abs(fft(dechirp_1_2));
dechirp_2_fft = abs(fft(dechirp_2));
dechirp_3_fft = abs(fft(dechirp_3));
dechirp_1_1_fft = dechirp_1_1_fft(1:fft_x) + dechirp_1_1_fft(dine-fft_x+1:dine);
dechirp_1_2_fft = dechirp_1_2_fft(1:fft_x) + dechirp_1_2_fft(dine-fft_x+1:dine);
dechirp_2_fft = dechirp_2_fft(1:fft_x) + dechirp_2_fft(dine-fft_x+1:dine);
dechirp_3_fft = dechirp_3_fft(1:fft_x) + dechirp_3_fft(dine-fft_x+1:dine);
subplot(4,1,1);
plot(dechirp_1_1_fft);
subplot(4,1,2);
plot(dechirp_1_2_fft);
subplot(4,1,3);
plot(dechirp_2_fft);
subplot(4,1,4);
plot(dechirp_3_fft);
% diff_1 = sum(abs(dechirp_1_fft - dechirp_2_fft));
% diff_2 = sum(abs(dechirp_1_fft - dechirp_3_fft));
% [~, max_bin_1] = max(dechirp_2_fft);
% [~, max_bin_2] = max(dechirp_3_fft);
% fre_dif = (max_bin_1 - max_bin_2)*(500000/200);

toc;
fclose all;
