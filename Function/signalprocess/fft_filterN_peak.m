function [pos_array, peak_array] = fft_filterN_peak(lora_set, fft_result, N)
    len = size(fft_result, 2);
    [peak,pos] = sort(fft_result, 2, 'descend');         % 对FFT进行排序
    pos_array = zeros(1, N);
    peak_array = zeros(1, N);
    pos_array(1) = pos(1);
    peak_array(1) = peak(1);
    count = 1;
    for i = 2:len
        temp1 = pos(i);
        for k = 1:count
            temp2 = pos(k);
            dif = abs(temp2 - temp1);
            if dif < len*lora_set.leakage_width1 || dif > len*lora_set.leakage_width2
                break;
            end
            if k == count 
                count = count + 1;
                pos_array(count) = temp1;
                peak_array(count) = peak(i);
            end
        end
        if count == N
            break;
        end
    end