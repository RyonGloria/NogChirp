function [lora_set] = readLoraSet(json_name, sf, bw, samples_rate, Config_Path)

    % Check if Config_Path is provided, otherwise use default value
    if nargin < 5 || isempty(Config_Path)
        Config_Path = '.\Config\';                               % 设置配置文件所在的文件夹路径
    end
    % 读取配置信息
    Setting_File = dir(fullfile(Config_Path, json_name));           % 找到配置文件夹下符合名字的文件
    Setting_File_Path = strcat(Config_Path, Setting_File.name);

    Setting_file = fopen(Setting_File_Path,'r');
    setting = jsondecode(fscanf(Setting_file,'%s'));                % 解析json格式变量
    if exist('sf', 'var')
        lora_set.sf = sf;
    else
        lora_set.sf = setting.captures.lora_sf;                         % 设置接收数据包的lora_set.sf
    end
    if exist('bw', 'var')
        lora_set.bw = bw;
    else
        lora_set.bw = setting.captures.lora_bw;                         % 设置接收数据包的lora_set.sf
    end
%     lora_set.sample_rate = setting.global.core_sample_rate;         % 设置接收数据包的采样率
%     lora_set.sample_rate = 2e6;
    if(~exist('samples_rate','var'))
        lora_set.sample_rate = setting.global.core_sample_rate;  % 如果未出现该变量，则对其进行赋值
    else
        lora_set.sample_rate = samples_rate;
    end

    lora_set.Pkg_length = setting.captures.lora_pkg_length;         % 设置接收数据包的长度
    lora_set.dine = lora_set.sample_rate*bitshift(1,lora_set.sf)/lora_set.bw;    % 根据lora_set.sf和lora_set.bw计算出一个chirp包含的采样点个数
    lora_set.fft_x = 2^lora_set.sf;                               % 根据lora_set.dine计算出包含lora_set.bw所需的FFT点数
    % pass_arg为matlab fir1滤波器中的参数，根据带宽调整滤波器的截止频率
    if lora_set.bw == 125e3
        lora_set.pass_arg = 0.05;
    elseif lora_set.bw == 250e3
        lora_set.pass_arg = 0.075;
    elseif lora_set.bw == 500e3
        lora_set.pass_arg = 0.15;
    end
    % preamble的数目
    lora_set.Preamble_length = 8;
    % 用于过滤FFT中旁瓣影响的参数，具体为fft_x * leakage_width1的范围内只允许存在一个峰，在这个范围内其他峰会被忽略掉
%     leakage_width_array = [0.05,0.01,0.015,0.01,0.01,0.01];
%     leakage_width_array = [0.05,0.01,0.01,0.005,0.01,0.01];
    leakage_width_array = [0.05,0.01,0.01,0.005,0.01,0.01];
    % 需要过滤出的峰的数目，2表示一个FFT需要过滤出的峰的数目为2+1个
    lora_set.filter_num = 2;
    if lora_set.sf > 6
        lora_set.leakage_width1 = leakage_width_array(lora_set.sf-6);
        lora_set.leakage_width2 = 1-lora_set.leakage_width1;
    end
    % 信道数目
    lora_set.channel = 4;
    % 选择的信道
    lora_set.channel_choice = [0,1,2,3];
    % 选择的信道的数目
    lora_set.channel_choice_num = 4;
    % preanble所在的信道（未使用）
    lora_set.Preamble_channel = 1;
    % 表示一个fft点对应多少个采样点
    lora_set.factor = lora_set.dine/lora_set.fft_x;

    fclose all;