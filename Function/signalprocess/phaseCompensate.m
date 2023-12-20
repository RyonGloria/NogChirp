function [dechirp_fft] = phaseCompensate(dechirp, loraSet)
    fftX = loraSet.fft_x;
    dine = loraSet.dine;
    dechirp_fft = fft(dechirp, dine);
    max_peak = 0;
    for pahse = 0:2*pi/180*128:2*pi
        dechirp_fft_result = dechirp_fft(1:fftX).*exp(1i*pahse) + dechirp_fft(dine-fftX+1:dine);
        dechirp_fft_result_abs = abs(dechirp_fft_result);
        max_value = max(dechirp_fft_result_abs(85:90));
        if max_value > max_peak
            max_peak = max_value;
            dechirp_fft_tmp = dechirp_fft_result_abs;
        end
    end
    dechirp_fft = dechirp_fft_tmp;