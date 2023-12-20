% 画出所选特定文件的相位图
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

Config_Path = '.\config\';                                       % 设置配置文件所在路径
% Verification_path = strcat(Config_Path,'node1_sf10.txt');        % bin值验证文件
Setting_File = dir(fullfile(Config_Path,'sf10_BW500.json'));     % 配置文件
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

fid_1=fopen('D:\CHchirp_samples\align_windows\USRP\samples_2.sigmf-data','rb');
[A_1]=fread(fid_1,'float32')';
A_1_length = size(A_1,2);
G0 = A_1(1:2:A_1_length-1) + A_1(2:2:A_1_length)*1i;            % 将float数组转换成复数数组
[d_downchirp, d_upchirp] = build_idealchirp(lora_set);
% [cfo, winoffset] = nscale_get_cfo(G0, lora_set);
[cfo, windows_offset] = get_cfo_winoff(G0, lora_set, d_downchirp, d_upchirp, 2);
G0 = circshift(G0, -round(windows_offset));
[d_downchirp_cfo, d_upchirp_cfo] = rebuild_idealchirp_cfo(lora_set, cfo);

dine = lora_set.dine;
fft_x = lora_set.fft_x;
% G0_tmp = G0(20.75*dine+1 : 21.75*dine);
% 
% dechirp = G0_tmp .* d_downchirp_cfo;
% dechirp_fft = fft(dechirp, lora_set.dine);

% result = zeros(0);
% result_bin = zeros(0);
% for ang = 0 : 2*pi/360 : 2*pi-2*pi/360
%     dechirp_fft_merge = dechirp_fft(1:fft_x)*exp(ang*1j) + dechirp_fft(fft_x+1:dine);
%     [tmp, tmp_pos] = max(dechirp_fft_merge);
%     result = [result, tmp];
%     result_bin = [result_bin, tmp_pos];
% end
% [~, tmp_pos] = max(result);
% if tmp_pos < 60 || tmp_pos > 120
%     G0 = circshift(G0, 1);
% end

% dechirp
i = 20;
G0_tmp_1 = G0(i*dine+0.25*dine+1 : i*dine+1.25*dine);
G0_tmp_2 = G0(i*dine+1.25*dine+1 : i*dine+2.25*dine);
dechirp = [G0_tmp_1.*d_downchirp_cfo, G0_tmp_2.*d_downchirp_cfo];

% 计算相位差
angle_tmp = zeros(1,dine*2);
for count = 1:dine*2-1
   iphase_1 = angle(dechirp(count));
   iphase_2 = angle(dechirp(count+1));
   while (iphase_2 - iphase_1) >  pi
       iphase_2 = iphase_2 - 2*pi;
   end
   while (iphase_2 - iphase_1) <  -pi
       iphase_2 = iphase_2 + 2*pi;
   end
   angle_tmp(count) = iphase_2 - iphase_1;
end
angle_tmp(dine*2) = angle_tmp(dine*2-1);
% 画出相位图
plot(angle_tmp,'');

% 计算chirp bin推算相位跳变的位置，计算相位均值和极差
% dechirp_1 = dechirp(1:dine);
% dechirp_2 = dechirp(dine+1:end);
% 
% dechirp_1_fft = abs(fft(dechirp_1));
% dechirp_2_fft = abs(fft(dechirp_2));
% dechirp_1_fft = dechirp_1_fft(1:fft_x) + dechirp_1_fft(dine-fft_x+1:dine);
% dechirp_2_fft = dechirp_2_fft(1:fft_x) + dechirp_2_fft(dine-fft_x+1:dine);
% [~, chirp1_bin] = max(dechirp_1_fft);
% [~, chirp2_bin] = max(dechirp_2_fft);
% chirp1_phasejump_pos = (fft_x-chirp1_bin)*2;
% chirp2_phasejump_pos = (fft_x-chirp2_bin)*2+dine;
% if chirp1_phasejump_pos < dine-20 || chirp2_phasejump_pos > dine+30
% %     chirp1_left_phase = mean(angle_tmp(30:chirp1_phasejump_pos-20));
%     chirp1_right_phase = mean(angle_tmp(chirp1_phasejump_pos+30:dine-20));
%     chirp2_left_phase = mean(angle_tmp(dine+30:chirp2_phasejump_pos-20));
% %     chirp2_right_phase = mean(angle_tmp(chirp2_phasejump_pos+30:dine*2-20));
%     chirp_minmax_lr_phase = zeros(1,8);
% %     chirp_minmax_lr_phase(1) = max(angle_tmp(30:chirp1_phasejump_pos-20));
% %     chirp_minmax_lr_phase(2) = min(angle_tmp(30:chirp1_phasejump_pos-20));
%     chirp_minmax_lr_phase(3) = max(angle_tmp(chirp1_phasejump_pos+30:dine-20));
%     chirp_minmax_lr_phase(4) = min(angle_tmp(chirp1_phasejump_pos+30:dine-20));
%     chirp_minmax_lr_phase(5) = max(angle_tmp(dine+30:chirp2_phasejump_pos-20));
%     chirp_minmax_lr_phase(6) = min(angle_tmp(dine+30:chirp2_phasejump_pos-20));
% %     chirp_minmax_lr_phase(7) = max(angle_tmp(chirp2_phasejump_pos+30:dine*2-20));
% %     chirp_minmax_lr_phase(8) = min(angle_tmp(chirp2_phasejump_pos+30:dine*2-20));
%     tmp = (angle_tmp(dine-20:dine+30) > chirp_minmax_lr_phase(3)) | (angle_tmp(dine-20:dine+30) < chirp_minmax_lr_phase(4));
%     tmp2 = (angle_tmp(dine-20:dine+30) > chirp_minmax_lr_phase(5)) | (angle_tmp(dine-20:dine+30) < chirp_minmax_lr_phase(6));
%     tmp3 = tmp & tmp2;
%     sum(tmp3)
% end