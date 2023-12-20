function autocorr = detect_preamble_autocorr(samples, lora_set)
    % samples: 输入序列
    window = lora_set.dine;
    chirp1 = samples(1:window);
    chirp2 = samples(window+1: 2*window);
    magsq_chirp1 = abs(chirp1) .^ 2;
    magsq_chirp2 = abs(chirp2) .^ 2;
    energy_chirp1 = sum(magsq_chirp1);
    energy_chirp2 = sum(magsq_chirp2);
    dot_product = conj(chirp1) .* chirp2;
    dot_product = sum(dot_product);
    
    autocorr = abs(dot_product / (sqrt(energy_chirp1 * energy_chirp2)));
end