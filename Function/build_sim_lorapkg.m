% 根据参数生成仿真lora信号
function [sim_samples] = build_sim_lorapkg(lora_set, payload_bin)
    % 设置参数
    sf = lora_set.sf;
    N = 2^sf;
    bw = lora_set.bw;
    dine = lora_set.dine;
    samples_rate = lora_set.sample_rate;
    os_factor = samples_rate / bw;
    C = lora_set.channel;   % channel number
    Cs = lora_set.channel_choice;
    pre_ch = lora_set.Preamble_channel;
    Cs_num = lora_set.channel_choice_num;
    
    bin = 0;
    n_fold = N*os_factor - bin*os_factor;
    % build upchirp matrix method
    n_org = 0:1:n_fold-1;
    n_jump = n_fold:1:N*os_factor-1;
    upchirp_org = (1+0*1i)*exp(2*pi*1i.*(n_org.*n_org./(2*N)/power(os_factor,2) + (2*pre_ch-C)*0.5.*n_org./os_factor + (bin/N).*n_org./os_factor));
    upchirp_jump = (1+0*1i)*exp(2*pi*1i.*(n_jump.*n_jump./(2*N)/power(os_factor,2) + (2*pre_ch-C)*0.5.*n_jump./os_factor + (bin/N-1).*n_jump./os_factor));
    upchirp_pre = [upchirp_org, upchirp_jump];

    % build downchirp matrix method
    pre_ch = C - pre_ch - 1;
    n_org = 0:1:n_fold-1;
    n_jump = n_fold:1:N*os_factor-1;
    upchirp_org = (1+0*1i)*exp(2*pi*1i.*(n_org.*n_org./(2*N)/power(os_factor,2) + (2*pre_ch-C)*0.5.*n_org./os_factor + (bin/N).*n_org./os_factor));
    upchirp_jump = (1+0*1i)*exp(2*pi*1i.*(n_jump.*n_jump./(2*N)/power(os_factor,2) + (2*pre_ch-C)*0.5.*n_jump./os_factor + (bin/N-1).*n_jump./os_factor));
    upchirp_tmp = [upchirp_org, upchirp_jump];
    downchirp_sfd = conj(upchirp_tmp);
    
    % build sync_symbol matrix method
    sync_symbol = zeros(1, N*os_factor*2);
    for bin_tmp = 1:2
        pre_ch = lora_set.Preamble_channel;
        bin = bin_tmp*8;
        n_fold = N*os_factor - bin*os_factor;

        n_org = 0:1:n_fold-1;
        n_jump = n_fold:1:N*os_factor-1;
        upchirp_org = (1+0*1i)*exp(2*pi*1i.*(n_org.*n_org./(2*N)/power(os_factor,2) + (2*pre_ch-C)*0.5.*n_org./os_factor + (bin/N).*n_org./os_factor));
        upchirp_jump = (1+0*1i)*exp(2*pi*1i.*(n_jump.*n_jump./(2*N)/power(os_factor,2) + (2*pre_ch-C)*0.5.*n_jump./os_factor + (bin/N-1).*n_jump./os_factor));
        sync_symbol((bin_tmp-1)*N*os_factor+1:bin_tmp*N*os_factor) = [upchirp_org, upchirp_jump];
    end
    pre_SFD = [repmat(upchirp_pre,1,8), sync_symbol, repmat(downchirp_sfd,1,2), downchirp_sfd(1:0.25*dine)];
    
    % build payload_symbol with payload bin by matrix method
    Cs_cal = zeros(1, N*os_factor);  % channel choice caculate array
    for i = 1:Cs_num
        Cs_cal((i-1)/Cs_num*N*os_factor+1:i/Cs_num*N*os_factor) = Cs(i);
    end
    pay_symbol = zeros(1, dine*length(payload_bin));
    for bin_index = 1:length(payload_bin)
        bin = payload_bin(bin_index);
        n_fold = N*os_factor - bin*os_factor;

        n_org = 0:1:n_fold-1;
        n_jump = n_fold:1:N*os_factor-1;
        upchirp_org = (1+0*1i)*exp(2*pi*1i.*(n_org.*n_org./(2*N)/power(os_factor,2) + (2*Cs_cal(n_org+1)-C)*0.5.*n_org./os_factor + (bin/N).*n_org./os_factor));
        upchirp_jump = (1+0*1i)*exp(2*pi*1i.*(n_jump.*n_jump./(2*N)/power(os_factor,2) + (2*Cs_cal(n_jump+1)-C)*0.5.*n_jump./os_factor + (bin/N-1).*n_jump./os_factor));
        pay_symbol((bin_index-1)*dine+1:bin_index*dine) = [upchirp_org, upchirp_jump];
    end
    sim_samples = [pre_SFD, pay_symbol];
