% 根据cfo生成解调payload的downchirp
function [d_downchirp_cfo] = rebuild_idealchirp_cfo_pay(lora_set, cfo)
    sf = lora_set.sf;
    N = 2^sf;
    bw = lora_set.bw;
    samples_rate = lora_set.sample_rate;
    os_factor = samples_rate / bw;
    
    C = lora_set.channel;   % channel number
    Cs = lora_set.channel_choice;
    Cs_num = lora_set.channel_choice_num;

    Cs_cal = zeros(1, N*os_factor);  % channel choice caculate array
    for i = 1:Cs_num
        Cs_cal((i-1)/Cs_num*N*os_factor+1:i/Cs_num*N*os_factor) = Cs(i);
    end
    bin = 0;
    n_fold = N*os_factor - bin*os_factor;

    % build downchirp matrix method
    n_org = 0:1:n_fold-1;
    n_jump = n_fold:1:N*os_factor-1;
%     upchirp_org = (1+0*1i)*exp(2*pi*1i.*(n_org.*n_org./(2*N)/power(os_factor,2) + (6*Cs_cal(n_org+1)-3*C+1-cfo*4/bw)*0.25.*n_org./os_factor + (bin/N).*n_org./os_factor));
%     upchirp_jump = (1+0*1i)*exp(2*pi*1i.*(n_jump.*n_jump./(2*N)/power(os_factor,2) + (6*Cs_cal(n_jump+1)-3*C+1-cfo*4/bw)*0.25.*n_jump./os_factor + (bin/N).*n_jump./os_factor));
    upchirp_org = (1+0*1i)*exp(2*pi*1i.*(n_org.*n_org./(2*N)/power(os_factor,2) + (2*Cs_cal(n_org+1)-C-cfo*2/bw)*0.5.*n_org./os_factor + (bin/N).*n_org./os_factor));
    upchirp_jump = (1+0*1i)*exp(2*pi*1i.*(n_jump.*n_jump./(2*N)/power(os_factor,2) + (2*Cs_cal(n_jump+1)-C-cfo*2/bw)*0.5.*n_jump./os_factor + (bin/N).*n_jump./os_factor));
    upchirp = [upchirp_org, upchirp_jump];
    d_downchirp_cfo = conj(upchirp);

    % build downchirp loop method, same as gr-lora_sdr
    % bin = 0;
    % n_fold = N*os_factor - bin*os_factor;
    % upchirp_tmp = zeros(1, N*os_factor);
    % for n = 0:1:N*os_factor-1
    %     Cs_choice = floor(n*Cs_num/(N*os_factor));
    %     if(n<n_fold)
    %         upchirp_tmp(n+1) = (1+0*1i)*exp(2*pi*1i*(n*n/(2*N)/power(os_factor,2) + (6*Cs(Cs_choice+1)-3*C+1-cfo*4/bw)*0.25*n/os_factor + (bin/N)*n/os_factor));
    %     else 
    %         upchirp_tmp(n+1) = (1+0*1i)*exp(2*pi*1i*(n*n/(2*N)/power(os_factor,2) + (6*Cs(Cs_choice+1)-3*C+1-cfo*4/bw)*0.25*n/os_factor + (bin/N-1)*n/os_factor));
    %     end
    % end
    % d_downchirp_cfo = conj(upchirp_tmp);

    