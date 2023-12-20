function [d_upchirp, d_downchirp] = build_idealchirp_tmp(lora_set, bin)
    sf = lora_set.sf;
    N = 2^sf;
    bw = lora_set.bw;
    samples_rate = lora_set.sample_rate;
    os_factor = samples_rate / bw;
    
    C = lora_set.channel;   % channel number
    Cs = lora_set.channel_choice;
    Cs_num = lora_set.channel_choice_num;

    n_fold = N*os_factor - bin*os_factor;
    d_upchirp = zeros(1, N*os_factor);
    for n = 0:1:N*os_factor-1
        Cs_choice = floor(n*Cs_num/(N*os_factor));
        if(n<n_fold)
            d_upchirp(n+1) = (1+0*1i)*exp(2*pi*1i*(n*n/(2*N)/power(os_factor,2) + (2*Cs(Cs_choice+1)-C)*0.5*n/os_factor + (bin/N)*n/os_factor));
        else 
            if Cs(Cs_choice+1) == C-1
                d_upchirp(n+1) = (1+0*1i)*exp(2*pi*1i*(n*n/(2*N)/power(os_factor,2) + (2*Cs(Cs_choice+1)-C)*0.5*n/os_factor + (bin/N-C)*n/os_factor));
            else
                d_upchirp(n+1) = (1+0*1i)*exp(2*pi*1i*(n*n/(2*N)/power(os_factor,2) + (2*Cs(Cs_choice+1)-C)*0.5*n/os_factor + (bin/N)*n/os_factor));
            end
        end
    end
    d_downchirp = conj(d_upchirp);