addpath('..');
rf_freq = 470e6;    % carrier frequency, used to correct clock drift
sf = 7;             % spreading factor
bw = 125e3;         % bandwidth
fs = 1e6;           % sampling rate

phy = LoRaPHY(rf_freq, sf, bw, fs);
phy.has_header = 1;         % explicit header mode
phy.cr = 4;                 % code rate = 4/8 (1:4/5 2:4/6 3:4/7 4:4/8)
phy.crc = 1;                % enable payload CRC checksum
phy.preamble_len = 8;       % preamble: 8 basic upchirps

% Encode payload [1 2 3 4 5]
% symbols = phy.encode((1:5)');
% fprintf("[encode] symbols:\n");
% disp(symbols);
symbols = (1:10:10*10)';
% Baseband Modulation
sig = phy.modulate(symbols);
% sig1 = [sig; zeros(phy.sample_num*4*15.5,1)];
% sig2 = [zeros(phy.sample_num*4*15.5,1); sig];
% save the complex signal to file
% samples = [sig1 + 0.3*sig2; zeros(phy.sample_num*4*40,1)];
% LoRaPHY.write(samples, './signals');
sig1 = [sig; zeros(phy.sample_num*4*40,1)];
sig2 = [zeros(phy.sample_num*4*16.25,1); sig; zeros(phy.sample_num*4*23.75,1)];
samples = sig1 + 0.3*sig2;
LoRaPHY.write(samples, 'E:\Pyramid_samples\conflict\conflict.sigmf-data');
