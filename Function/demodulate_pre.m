% 解调得到输入流的preamble的bin值
function [max_pos_array] = demodulate_pre(G0, lora_set, upchirp, downchirp)
    dine = lora_set.dine;
    fft_x = lora_set.fft_x;
    tmp = [repmat(downchirp, 1, lora_set.Preamble_length+2), repmat(upchirp, 1, 2)];
    len_all = length(tmp) / lora_set.dine;   % 计算需要解调的窗口数目
    tmp = reshape(tmp, [dine, len_all]).';   % 重组chirp的形式，一行为一个chirp，列数为chirp的数目
    samples = reshape(G0(1:len_all*dine), [dine, len_all]).';

    samples_dechirp = samples .* tmp;
    samples_fft = abs(fft(samples_dechirp,dine,2));
    samples_fft_merge = samples_fft(:,1:fft_x) + samples_fft(:,dine-fft_x+1:dine);
    [~, max_pos_array] = max(samples_fft_merge, [], 2);