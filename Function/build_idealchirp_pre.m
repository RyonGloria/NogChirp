% 生成解调preamble的idealchirp
function [downchirp, upchirp] = build_idealchirp_pre(lora_set)
    sf = lora_set.sf;
    N = 2^sf;
    bw = lora_set.bw;
    samples_rate = lora_set.sample_rate;
    os_factor = samples_rate / bw;
    
    C = lora_set.channel;   % channel number
    pre_ch = lora_set.Preamble_channel;
    pre_ch = C - pre_ch - 1;

    bin = 0;
    n_fold = N*os_factor - bin*os_factor;

    % build upchirp matrix method
    n_org = 0:1:n_fold-1;
    n_jump = n_fold:1:N*os_factor-1;
%     upchirp_org = (1+0*1i)*exp(2*pi*1i.*(n_org.*n_org./(2*N)/power(os_factor,2) + (6*pre_ch-3*C+1)*0.25.*n_org./os_factor + (bin/N).*n_org./os_factor));
%     upchirp_jump = (1+0*1i)*exp(2*pi*1i.*(n_jump.*n_jump./(2*N)/power(os_factor,2) + (6*pre_ch-3*C+1)*0.25.*n_jump./os_factor + (bin/N).*n_jump./os_factor));
    upchirp_org = (1+0*1i)*exp(2*pi*1i.*(n_org.*n_org./(2*N)/power(os_factor,2) + (2*pre_ch-C)*0.5.*n_org./os_factor + (bin/N).*n_org./os_factor));
    upchirp_jump = (1+0*1i)*exp(2*pi*1i.*(n_jump.*n_jump./(2*N)/power(os_factor,2) + (2*pre_ch-C)*0.5.*n_jump./os_factor + (bin/N).*n_jump./os_factor));
    upchirp = [upchirp_org, upchirp_jump];

    % build downchirp matrix method
    pre_ch = lora_set.Preamble_channel;
    n_org = 0:1:n_fold-1;
    n_jump = n_fold:1:N*os_factor-1;
%     upchirp_org = (1+0*1i)*exp(2*pi*1i.*(n_org.*n_org./(2*N)/power(os_factor,2) + (6*pre_ch-3*C+1)*0.25.*n_org./os_factor + (bin/N).*n_org./os_factor));
%     upchirp_jump = (1+0*1i)*exp(2*pi*1i.*(n_jump.*n_jump./(2*N)/power(os_factor,2) + (6*pre_ch-3*C+1)*0.25.*n_jump./os_factor + (bin/N).*n_jump./os_factor));
    upchirp_org = (1+0*1i)*exp(2*pi*1i.*(n_org.*n_org./(2*N)/power(os_factor,2) + (2*pre_ch-C)*0.5.*n_org./os_factor + (bin/N).*n_org./os_factor));
    upchirp_jump = (1+0*1i)*exp(2*pi*1i.*(n_jump.*n_jump./(2*N)/power(os_factor,2) + (2*pre_ch-C)*0.5.*n_jump./os_factor + (bin/N).*n_jump./os_factor));
    upchirp_tmp = [upchirp_org, upchirp_jump];
    downchirp = conj(upchirp_tmp);
    
    % build upchirp loop method, same as gr-lora_sdr
    % bin = 0;
    % n_fold = N*os_factor - bin*os_factor;
    % upchirp = zeros(1, N*os_factor);
    % for n = 0:1:N*os_factor-1
    %     if(n<n_fold)
    %         upchirp(n+1) = (1+0*1i)*exp(2*pi*1i*(n*n/(2*N)/power(os_factor,2) + (6*pre_ch-3*C+1)*0.25*n/os_factor + (bin/N)*n/os_factor));
    %     else 
    %         upchirp(n+1) = (1+0*1i)*exp(2*pi*1i*(n*n/(2*N)/power(os_factor,2) + (6*pre_ch-3*C+1)*0.25*n/os_factor + (bin/N-1)*n/os_factor));
    %     end
    % end

    % % build downchirp loop method, same as gr-lora_sdr
    % pre_ch = lora_set.Preamble_channel;
    % bin = 0;
    % n_fold = N*os_factor - bin*os_factor;
    % upchirp_tmp = zeros(1, N*os_factor);
    % for n = 0:1:N*os_factor-1
    %     if(n<n_fold)
    %         upchirp_tmp(n+1) = (1+0*1i)*exp(2*pi*1i*(n*n/(2*N)/power(os_factor,2) + (6*pre_ch-3*C+1)*0.25*n/os_factor + (bin/N)*n/os_factor));
    %     else 
    %         upchirp_tmp(n+1) = (1+0*1i)*exp(2*pi*1i*(n*n/(2*N)/power(os_factor,2) + (6*pre_ch-3*C+1)*0.25*n/os_factor + (bin/N-1)*n/os_factor));
    %     end
    % end
    % downchirp = conj(upchirp_tmp);