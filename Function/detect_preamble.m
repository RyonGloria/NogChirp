% 检测输入流中是否存在preamble，返回检测到的preamble窗口位置
function [Preamble_start_pos] = detect_preamble(G0, lora_set, idealchirp, detect_length)
    % 检测Preamble
    % 获取Preamble，dechirp，计算FFT并合并FFT
    dine = lora_set.dine;
    fft_x = lora_set.fft_x;
    Preamble_length = lora_set.Preamble_length;
    if ~exist("detect_length")
        detect_length = Preamble_length*2;
    end
    samples = reshape(G0(1:detect_length*dine), [dine, detect_length]).';
    samples_dechirp = samples .* idealchirp;
    samples_fft = abs(fft(samples_dechirp, dine, 2));
    samples_fft_merge = samples_fft(:,1:fft_x) + samples_fft(:,dine-fft_x+1:dine);
%     plot(samples_fft_merge(1,:));
%     [pos_array, peak_array] = fft_filterN_peak(lora_set, samples_fft_merge(1,:), 3);
    
    [fft_amp_max, ~] = max(samples_fft_merge,[],2);    % 获得FFT中的主峰的峰值强度和位置
    % 找到Preamble开始位置，如果不等于1，要将窗口位移
    Preamble_start_pos = find((fft_amp_max > 3*mean(samples_fft_merge,2)) == 1,1,'first');
    if Preamble_start_pos == []
        fprintf("Can't detect preamble!\n");
        Preamble_start_pos = 999;
        return;
    end

