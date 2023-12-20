function [bin_array, channel_choice_array] = get_channel_bin(G0, lora_set, cfo, payload_length)
    dine = lora_set.dine;
    fft_x = lora_set.fft_x;
    Ch_num = lora_set.channel_choice_num;
    start_pos = lora_set.Preamble_length+2+2.25;
    G0_pay = G0(start_pos*dine+1:(start_pos+payload_length)*dine);
    % 生成所有信道的downchirp
    downchirp_array = zeros(Ch_num, dine);
    for i = 1:Ch_num
        lora_set.channel_choice = repmat(i-1, 1, Ch_num);
        downchirp_array(i, :) = rebuild_idealchirp_cfo_pay(lora_set, cfo);
    end
    % 对每个chirp遍历所有信道可能，找到正确的bin和信道策略
    bin_array = zeros(1, payload_length);
    channel_choice_array = zeros(payload_length, Ch_num);
    for i = 1:payload_length
        samples = G0_pay((i-1)*dine+1:i*dine);
        % 获得bin值
        lora_set.channel_choice = zeros(1,Ch_num);
        downchirp_tmp = rebuild_idealchirp_cfo_pay(lora_set, cfo);
        dechirp = samples.*downchirp_tmp;
        [~,bin_tmp] = max(abs(fft(dechirp)));
        bin = mod(bin_tmp, fft_x);
        if bin == 0
            bin = 1;
        end
        bin_array(i) = bin;
        % 根据bin判断chirp跳转的位置
        chirp_jump_pos = Ch_num - floor(bin/(fft_x/Ch_num));
        sub_length = dine / Ch_num;
        choice_all = zeros(Ch_num, Ch_num);
        for chirp_sub = 1:Ch_num
            samples_sub = samples((chirp_sub-1)*sub_length+1 : chirp_sub*sub_length);
            for channel_sub = 1:Ch_num
                downchirp_sub = downchirp_array(channel_sub, (chirp_sub-1)*sub_length+1 : chirp_sub*sub_length);
                dechirp_sub = samples_sub.*downchirp_sub;
                fft_sub = abs(fft(dechirp_sub, dine));
                if chirp_sub == chirp_jump_pos
                    bin_mod = mod(bin-fft_x, dine);
                    choice_all(chirp_sub, channel_sub) = fft_sub(bin)+fft_sub(bin_mod);
                elseif chirp_sub > chirp_jump_pos
                    bin_mod = mod(bin-fft_x, dine);
                    choice_all(chirp_sub, channel_sub) = fft_sub(bin_mod);
                else
                    choice_all(chirp_sub, channel_sub)= fft_sub(bin);
                end
            end
        end
        [~, tmp] = max(choice_all,[],2);
        channel_choice_array(i,:) = tmp-1;
    end