% 补零计算lora的cfo和time_offset
function [cfo, windows_offset] = get_cfo_winoff(G0, lora_set, downchirp, upchirp, Preamble_num, decim, ifCHchip)
    % 计算主峰的CFO(需要补零操作)
    % 对Preamble阶段的FFT峰值进行排序，得到前filter的峰
    zeropadding_size = decim;                   % 设置补零的数量，这里的decim表示，补上decim-1倍窗口的零，计算FFT时一共是decim倍的窗口（decim+1）
    d_sf = lora_set.sf;
    d_bw = lora_set.bw;
    dine = lora_set.dine;
    fft_x = lora_set.fft_x;

    dine_zeropadding = dine*zeropadding_size*2^(10-d_sf);
    fft_x_zeropadding = fft_x*zeropadding_size*2^(10-d_sf);

    filter_num = lora_set.filter_num;
%     Preamble_length = lora_set.Preamble_length;
    leakage_width1 = lora_set.leakage_width1;
    leakage_width2 = lora_set.leakage_width2;
    samples = reshape(G0(1:Preamble_num*dine),[dine,Preamble_num]).';
    samples_fft = abs(fft(samples .* downchirp, dine_zeropadding, 2));
%     samples_fft_merge = [samples_fft(:,1:4096) + samples_fft(:,57345:61440) , samples_fft(:,61441:65536) + samples_fft(:,4097:8192)];
    samples_fft_merge = samples_fft(:, 1:fft_x_zeropadding) + samples_fft(:,dine_zeropadding-fft_x_zeropadding+1:dine_zeropadding);
    [peak,pos] = sort(samples_fft_merge(1:Preamble_num,:),2,'descend');         % 对FFT进行排序
    Peak_pos = zeros(size(pos,1),filter_num);
    Peak_amp = zeros(size(peak,1),filter_num);
    Peak_pos(:,1) = pos(:,1);
    Peak_amp(:,1) = peak(:,1);
    for row = 1:size(pos,1)
        temp_array = ones(1,size(pos,2));
        for list = 1:filter_num
            temp_array = temp_array & (abs(Peak_pos(row,list) - pos(row,:)) > fft_x_zeropadding*leakage_width1 & abs(Peak_pos(row,list) - pos(row,:)) < fft_x_zeropadding*leakage_width2);
            temp_num = find(temp_array==1,1,'first');
            Peak_pos(row,list+1) = pos(row,temp_num);
            Peak_amp(row,list+1) = peak(row,temp_num);
        end
    end
    
    %寻找与第一个窗口的峰（默认第一个窗口只有包1的峰）相近的峰，得到与其相近且重复次数最多的bin，记作Preamble的bin
    if Peak_pos(1) == Peak_pos(2)
        upchirp_ref = Peak_pos(1);
    else
        upchirp_ref = Peak_pos(2);
    end
    upchirp_index = abs(Peak_pos-upchirp_ref) < fft_x_zeropadding*leakage_width1 | abs(Peak_pos-upchirp_ref) > fft_x_zeropadding*leakage_width2;
    upchirp_bin = (Peak_pos(upchirp_index));
    upchirp_peak = mode(upchirp_bin);

    % 已知downchirp的位置，得到downchirp的bin
%     if ifCHchip == false
%         SFD_samples = reshape(G0((Preamble_length+2)*dine+1:(Preamble_length+4)*dine),[dine,2]).';
%     else
%         SFD_samples = reshape(G0((Preamble_length+2)*dine+1:(Preamble_length+4)*dine),[dine,2]).';
% %         SFD_samples = reshape(G0((Preamble_length+3)*dine+1:(Preamble_length+5)*dine),[dine,2]).';
%     end
    SFD_samples = G0((Preamble_num+2)*dine+1:(Preamble_num+3)*dine);
    SFD_samples_fft = abs(fft(SFD_samples .* upchirp, dine_zeropadding));
%     samples_fft_merge = [SFD_samples_fft(:,1:4096) + SFD_samples_fft(:,57345:61440) , SFD_samples_fft(:,61441:65536) + SFD_samples_fft(:,4097:8192)];
    samples_fft_merge = SFD_samples_fft(1:fft_x_zeropadding) + SFD_samples_fft(dine_zeropadding-fft_x_zeropadding+1:dine_zeropadding);
    [~, downchirp_peak] = max(samples_fft_merge);
%     downchirp_peak = SFD_pos_max;
    
    % 计算CFO和窗口偏移量
    if upchirp_peak + downchirp_peak < fft_x_zeropadding*0.5
        cfo_bin = upchirp_peak + downchirp_peak - 2;
        cfo = -cfo_bin/2/fft_x_zeropadding * d_bw;
        windows_offset = (downchirp_peak - upchirp_peak) / 2^(11-d_sf);
    elseif upchirp_peak + downchirp_peak > fft_x_zeropadding*1.5
        cfo_bin = upchirp_peak + downchirp_peak - fft_x_zeropadding*2 - 2;
        cfo = -cfo_bin/2/fft_x_zeropadding * d_bw;
        windows_offset = (downchirp_peak - upchirp_peak) / 2^(11-d_sf);
    else
        cfo_bin = upchirp_peak + downchirp_peak - fft_x_zeropadding - 2;
        cfo = -cfo_bin/2/fft_x_zeropadding * d_bw;
        windows_offset = (fft_x_zeropadding - (upchirp_peak - downchirp_peak)) / 2^(11-d_sf);
    end