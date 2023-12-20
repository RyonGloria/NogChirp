function [true_chirp, true_rate] = vertify_bin(subchirp_bin, true_bin)
    true_chirp = zeros(1, size(subchirp_bin, 1));
    subchirpNum = size(subchirp_bin, 2);
    for chirp = 1:size(subchirp_bin, 1)
        if sum(subchirp_bin(chirp, :) == ones(1, subchirpNum)*true_bin(chirp)) == subchirpNum
            true_chirp(chirp) = 1;
        end
    end
    true_rate = sum(true_chirp) / size(subchirp_bin, 1);