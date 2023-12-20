% 解调得到输入流的payload的bin值
function [max_pos_array] = demodulate_pay(G0, lora_set, downchirp, payload_length)
    dine = lora_set.dine;
    fft_x = lora_set.fft_x;
    start_pos = lora_set.Preamble_length+2+2.25;
    samples = reshape(G0(start_pos*dine+1:(start_pos+payload_length)*dine), [dine, payload_length]).';

    samples_dechirp = samples .* downchirp;
    samples_fft = abs(fft(samples_dechirp,dine,2));
    samples_fft_merge = samples_fft(:,1:fft_x) + samples_fft(:,dine-fft_x+1:dine);
    [~, max_pos_array] = max(samples_fft_merge, [], 2);