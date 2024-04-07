classdef LoRaPhyLay < handle & matlab.mixin.Copyable
    properties
        sf                        % spreading factor (7,8,9,10,11,12)
        bw                        % bandwidth (125kHz 250kHz 500kHz)
        fs                        % sampling frequency
        payload_len               % payload length
        cr                        % code rate: (1:4/5 2:4/6 3:4/7 4:4/8)
        has_header                % explicit header: 1, implicit header: 0
        crc                       % crc = 1 if CRC Check is enabled else 0
        ldr                       % ldr = 1 if Low Data Rate Optimization is enabled else 0
        crc_generator             % CRC generator with polynomial x^16+x^12+x^5+1
        header_checksum_matrix    % we use a 12 x 5 matrix to calculate header checksum
        whitening_seq             % whitening sequence
        hamming_decoding_en       % enable hamming decoding
        is_debug                  % set `true` for debug information
    end
    methods
        function obj = LoRaPhyLay(sf, bw, fs)
            obj.sf = sf;
            obj.bw = bw;
            obj.fs = fs;
            obj.has_header = 1;
            obj.crc = 1;
            obj.hamming_decoding_en = true;

            obj.whitening_seq = uint8([0xff, 0xfe, 0xfc, 0xf8, 0xf0, 0xe1, 0xc2, 0x85, 0xb, 0x17, 0x2f, 0x5e, 0xbc, 0x78, 0xf1, 0xe3, 0xc6, 0x8d, 0x1a, 0x34, 0x68, 0xd0, 0xa0, 0x40, 0x80, 0x1, 0x2, 0x4, 0x8, 0x11, 0x23, 0x47, 0x8e, 0x1c, 0x38, 0x71, 0xe2, 0xc4, 0x89, 0x12, 0x25, 0x4b, 0x97, 0x2e, 0x5c, 0xb8, 0x70, 0xe0, 0xc0, 0x81, 0x3, 0x6, 0xc, 0x19, 0x32, 0x64, 0xc9, 0x92, 0x24, 0x49, 0x93, 0x26, 0x4d, 0x9b, 0x37, 0x6e, 0xdc, 0xb9, 0x72, 0xe4, 0xc8, 0x90, 0x20, 0x41, 0x82, 0x5, 0xa, 0x15, 0x2b, 0x56, 0xad, 0x5b, 0xb6, 0x6d, 0xda, 0xb5, 0x6b, 0xd6, 0xac, 0x59, 0xb2, 0x65, 0xcb, 0x96, 0x2c, 0x58, 0xb0, 0x61, 0xc3, 0x87, 0xf, 0x1f, 0x3e, 0x7d, 0xfb, 0xf6, 0xed, 0xdb, 0xb7, 0x6f, 0xde, 0xbd, 0x7a, 0xf5, 0xeb, 0xd7, 0xae, 0x5d, 0xba, 0x74, 0xe8, 0xd1, 0xa2, 0x44, 0x88, 0x10, 0x21, 0x43, 0x86, 0xd, 0x1b, 0x36, 0x6c, 0xd8, 0xb1, 0x63, 0xc7, 0x8f, 0x1e, 0x3c, 0x79, 0xf3, 0xe7, 0xce, 0x9c, 0x39, 0x73, 0xe6, 0xcc, 0x98, 0x31, 0x62, 0xc5, 0x8b, 0x16, 0x2d, 0x5a, 0xb4, 0x69, 0xd2, 0xa4, 0x48, 0x91, 0x22, 0x45, 0x8a, 0x14, 0x29, 0x52, 0xa5, 0x4a, 0x95, 0x2a, 0x54, 0xa9, 0x53, 0xa7, 0x4e, 0x9d, 0x3b, 0x77, 0xee, 0xdd, 0xbb, 0x76, 0xec, 0xd9, 0xb3, 0x67, 0xcf, 0x9e, 0x3d, 0x7b, 0xf7, 0xef, 0xdf, 0xbf, 0x7e, 0xfd, 0xfa, 0xf4, 0xe9, 0xd3, 0xa6, 0x4c, 0x99, 0x33, 0x66, 0xcd, 0x9a, 0x35, 0x6a, 0xd4, 0xa8, 0x51, 0xa3, 0x46, 0x8c, 0x18, 0x30, 0x60, 0xc1, 0x83, 0x7, 0xe, 0x1d, 0x3a, 0x75, 0xea, 0xd5, 0xaa, 0x55, 0xab, 0x57, 0xaf, 0x5f, 0xbe, 0x7c, 0xf9, 0xf2, 0xe5, 0xca, 0x94, 0x28, 0x50, 0xa1, 0x42, 0x84, 0x9, 0x13, 0x27, 0x4f, 0x9f, 0x3f, 0x7f]');

            obj.header_checksum_matrix = gf([
                1 1 1 1 0 0 0 0 0 0 0 0
                1 0 0 0 1 1 1 0 0 0 0 1
                0 1 0 0 1 0 0 1 1 0 1 0
                0 0 1 0 0 1 0 1 0 1 1 1
                0 0 0 1 0 0 1 0 1 1 1 1
            ]);

            obj.crc_generator = comm.CRCGenerator('Polynomial','X^16 + X^12 + X^5 + 1');

            obj.init();
        end


        function init(obj)
            % Low Data Rate Optimization (LDRO) mode in LoRa
            % If the chirp peird is larger than 16ms, the least significant
            % two bits are considered unreliable and are neglected.
            if 2^(obj.sf)/obj.bw > 16e-3
                obj.ldr = 1;
            else
                obj.ldr = 0;
            end
        end

        function symbols = encode(obj, payload)
            % encode  Encode bytes to symbols
            %
            % input:
            %     payload: Payload of LoRa packet
            % output:
            %     symbols: A vector representing the symbols of the packet

            if obj.crc
                data = uint8([payload; obj.calc_crc(payload)]);
            else
                data = uint8(payload);
            end

            plen = length(payload);
            sym_num = obj.calc_sym_num(plen);
            % filling all symbols needs nibble_num nibbles
            nibble_num = obj.sf - 2 + (sym_num-8)/(obj.cr+4)*(obj.sf-2*obj.ldr);
            data_w = uint8([data; 255*ones(ceil((nibble_num-2*length(data))/2), 1)]);
            data_w(1:plen) = obj.whiten(data_w(1:plen));
            data_nibbles = uint8(zeros(nibble_num, 1));
            for i = 1:nibble_num
                idx = ceil(i/2);
                if mod(i, 2) == 1
                    data_nibbles(i) = bitand(data_w(idx), 0xf);
                else
                    data_nibbles(i) = bitshift(data_w(idx), -4);
                end
            end

            if obj.has_header
                header_nibbles = obj.gen_header(plen);
            else
                header_nibbles = [];
            end
            codewords = obj.hamming_encode([header_nibbles; data_nibbles]);

            % interleave
            % first 8 symbols use CR=4/8
            symbols_i = obj.diag_interleave(codewords(1:obj.sf-2), 8);
            ppm = obj.sf - 2*obj.ldr;
            rdd = obj.cr + 4;
            for i = obj.sf-1:ppm:length(codewords)-ppm+1
                symbols_i = [symbols_i; obj.diag_interleave(codewords(i:i+ppm-1), rdd)];
            end

            symbols = obj.gray_decoding(symbols_i);
        end

        function header_nibbles = gen_header(obj, plen)
            % gen_header  Generate a valid LoRa header
            %
            % input:
            %     plen: Payload length
            % output:
            %     header_nibbles: Header in nibbles

            header_nibbles = zeros(5, 1);
            header_nibbles(1) = bitshift(plen, -4);
            header_nibbles(2) = bitand(plen, 15);
            header_nibbles(3) = bitor(2*obj.cr, obj.crc);
            header_checksum = obj.header_checksum_matrix * gf(reshape(de2bi(header_nibbles(1:3), 4, 'left-msb')', [], 1));
            x = header_checksum.x;
            header_nibbles(4) = x(1);
            for i = 1:4
                header_nibbles(5) = bitor(header_nibbles(5), x(i+1)*2^(4-i));
            end
        end

        function checksum = calc_crc(obj, data)
            % calc_crc  Calculate payload CRC
            %
            % input:
            %     data: Data in bytes
            % output:
            %     checksum: CRC result

            switch length(data)
                case 0
                    checksum = [0; 0];
                case 1
                    checksum = [data(end); 0];
                case 2
                    checksum = [data(end); data(end-1)];
                otherwise
                    input = data(1:end-2);
                    seq = obj.crc_generator(reshape(logical(de2bi(input, 8, 'left-msb'))', [], 1));
                    checksum_b1 = bitxor(bi2de(seq(end-7:end)', 'left-msb'), data(end));
                    checksum_b2 = bitxor(bi2de(seq(end-15:end-8)', 'left-msb'), data(end-1));
                    checksum = [checksum_b1; checksum_b2];
            end
        end

        function data_w = whiten(obj, data)
            % whiten  Whitening process in LoRa
            %
            % input:
            %     data: Data in bytes
            % output:
            %     data_w: Data after whitening

            len = length(data);
            data_w = bitxor(data(1:len), obj.whitening_seq(1:len));
            obj.print_bin("Whiten", data_w);
        end

        % 汉明码编码
        function codewords = hamming_encode(obj, nibbles)
            % hamming_encode  Hamming encoding process in LoRa
            %
            % input:
            %     nibbles: Data in nibbles
            % output:
            %     codewords: Data after hamming encoding

            nibble_num = length(nibbles);
            codewords = uint8(zeros(nibble_num, 1));
            for i = 1:nibble_num
                nibble = nibbles(i);

                p1 = LoRaPHY.bit_reduce(@bitxor, nibble, [1 3 4]);
                p2 = LoRaPHY.bit_reduce(@bitxor, nibble, [1 2 4]);
                p3 = LoRaPHY.bit_reduce(@bitxor, nibble, [1 2 3]);
                p4 = LoRaPHY.bit_reduce(@bitxor, nibble, [1 2 3 4]);
                p5 = LoRaPHY.bit_reduce(@bitxor, nibble, [2 3 4]);
                if i <= obj.sf - 2
                    % the first SF-2 nibbles use CR=4/8
                    cr_now = 4;
                else
                    cr_now = obj.cr;
                end
                switch cr_now
                    case 1
                        codewords(i) = bitor(uint8(p4)*16, nibble);
                    case 2
                        codewords(i) = LoRaPHY.word_reduce(@bitor, [uint8(p5)*32 uint8(p3)*16 nibble]);
                    case 3
                        codewords(i) = LoRaPHY.word_reduce(@bitor, [uint8(p2)*64 uint8(p5)*32 uint8(p3)*16 nibble]);
                    case 4
                        codewords(i) = LoRaPHY.word_reduce(@bitor, [uint8(p1)*128 uint8(p2)*64 uint8(p5)*32 uint8(p3)*16 nibble]);
                    otherwise
                        % THIS CASE SHOULD NOT HAPPEN
                        error('Invalid Code Rate!');
                end
            end
        end

        function symbols_i = diag_interleave(obj, codewords, rdd)
            % diag_interleave  Diagonal interleaving
            %
            % input:
            %     codewords: Data in nibbles
            %     rdd: Bits with redundancy
            %          For example, code rate 4/5 means rdd = 5
            % output:
            %     symbols_i: Symbols after diagonal interleaving

            tmp = de2bi(codewords, rdd, 'right-msb');
            symbols_i = uint16(bi2de(cell2mat(arrayfun(@(x) circshift(tmp(:,x), 1-x), 1:rdd, 'un', 0))'));
            obj.print_bin("Interleave", symbols_i);
        end

        function symbols = gray_decoding(obj, symbols_i)
            % gray_decoding  Gray decoding
            %                `gray_decoding` is used in the ENCODING process
            %
            % input:
            %     symbols_i: Interleaved symbols
            % output:
            %     symbols: Final symbols to be modulated in a packet

            symbols = zeros(length(symbols_i), 1);
            for i = 1:length(symbols_i)
                num = uint16(symbols_i(i));
                mask = bitshift(num, -1);
                while mask ~= 0
                    num = bitxor(num, mask);
                    mask = bitshift(mask, -1);
                end
                if i <= 8 || obj.ldr
                    symbols(i) = mod(num * 4 + 1, 2^obj.sf);
                else
                    symbols(i) = mod(num + 1, 2^obj.sf);
                end
            end
        end

        function sym_num = calc_sym_num(obj, plen)
            % calc_sym_num  Calculate number of symbols
            %
            % input:
            %     plen: Payload length
            % output:
            %     sym_num: Number of symbols

            sym_num = double(8 + max((4+obj.cr)*ceil(double((2*plen-obj.sf+7+4*obj.crc-5*(1-obj.has_header)))/double(obj.sf-2*obj.ldr)), 0));
        end

        function [data_m, checksum_m] = decode(obj, symbols_m)
            % decode  Decode data from symbols
            %
            % input:
            %     symbols_m: A matrix of symbols to be decoded. Each column
            %                vector represents demodulated symbols from a
            %                LoRa packet.
            % output:
            %     data_m: A matrix of bytes representing the decoding
            %             result of `symbols_m`. The last two bytes are the
            %             decoded CRC checksum if CRC is enabled.
            %     checksum_m: A vector of checksum based on the decoded
            %                 payload.Checksum is empty if CRC is disabled.

            data_m = [];
            checksum_m = [];

            for pkt_num = 1:size(symbols_m, 2)
                % gray coding
                symbols_g = obj.gray_coding(symbols_m(:, pkt_num));

                % deinterleave
                codewords = obj.diag_deinterleave(symbols_g(1:8), obj.sf-2);
                if ~obj.has_header
                    nibbles = obj.hamming_decode(codewords, 8);
                else
                    % parse header
                    nibbles = obj.hamming_decode(codewords, 8);
                    obj.payload_len = double(nibbles(1)*16 + nibbles(2));
                    obj.crc = double(bitand(nibbles(3), 1));
                    obj.cr = double(bitshift(nibbles(3), -1));
                    % we only calculate header checksum on the first three nibbles
                    % the valid header checksum is considered to be 5 bits
                    % other 3 bits require further reverse engineering
                    header_checksum = [bitand(nibbles(4), 1); de2bi(nibbles(5), 4, 'left-msb')'];
                    header_checksum_calc = obj.header_checksum_matrix * gf(reshape(de2bi(nibbles(1:3), 4, 'left-msb')', [], 1));
                    if any(header_checksum ~= header_checksum_calc)
                        error('Invalid header checksum!');
                    end
                    nibbles = nibbles(6:end);
                end
                rdd = obj.cr + 4;
                for ii = 9:rdd:length(symbols_g)-rdd+1
                    codewords = obj.diag_deinterleave(symbols_g(ii:ii+rdd-1), obj.sf-2*obj.ldr);
                    % hamming decode
                    nibbles = [nibbles; obj.hamming_decode(codewords, rdd)];
                end

                % combine nibbles to bytes
                bytes = uint8(zeros(min(255, floor(length(nibbles)/2)), 1));
                for ii = 1:length(bytes)
                    bytes(ii) = bitor(uint8(nibbles(2*ii-1)), 16*uint8(nibbles(2*ii)));
                end

                % dewhitening
                len = obj.payload_len;
                if obj.crc
                    % The last 2 bytes are CRC16 checkcum
                    data = [obj.dewhiten(bytes(1:len)); bytes(len+1:len+2)];
                    % Calculate CRC checksum
                    checksum = obj.calc_crc(data(1:len));
                else
                    data = obj.dewhiten(bytes(1:len));
                    checksum = [];
                end
                data_m = [data_m data];
                checksum_m = [checksum_m checksum];
            end
        end

        function symbols = dynamic_compensation(obj, data)
            % dynamic_compensation  Compensate bin drift
            %
            % input:
            %     data: Symbols with bin drift
            % output:
            %     symbols: Symbols after bin calibration

            % compensate the bin drift caused by Sampling Frequency Offset (SFO)
            sfo_drift = (1 + (1:length(data))') * 2^obj.sf * obj.cfo / obj.rf_freq;
            symbols = mod(data - sfo_drift, 2^obj.sf);

            if obj.ldr
                bin_offset = 0;
                v_last = 1;

                for i = 1:length(symbols)
                    v = symbols(i);
                    bin_delta = mod(v-v_last, 4);
                    if bin_delta < 2
                        bin_offset = bin_offset - bin_delta;
                    else
                        bin_offset = bin_offset - bin_delta + 4;
                    end
                    v_last = v;
                    symbols(i) = mod(v+bin_offset, 2^obj.sf);
                end
            end
        end
        function symbols = gray_coding(obj, din)
            % gray_coding  Gray coding
            %              `gray_coding` is used in the DECODING process
            %
            % input:
            %     data: Symbols with bin drift
            % output:
            %     symbols: Symbols after bin calibration

            din(1:8) = floor(din(1:8)/4);
            if obj.ldr
                din(9:end) = floor(din(9:end)/4);
            else
                din(9:end) = mod(din(9:end)-1, 2^obj.sf);
            end
            s = uint16(din);
            symbols = bitxor(s, bitshift(s, -1));
            obj.print_bin("Gray Coding", symbols, obj.sf);
        end

        function codewords = diag_deinterleave(obj, symbols, ppm)
            % diag_deinterleave  Diagonal deinterleaving
            %
            % input:
            %     symbols: Symbols after gray coding
            %     ppm: Size with parity bits
            % output:
            %     codewords: Codewords after deinterleaving

            b = de2bi(symbols, double(ppm), 'left-msb');
            codewords = flipud(bi2de(cell2mat(arrayfun(@(x) ...
                circshift(b(x,:), [1 1-x]), (1:length(symbols))', 'un', 0))'));
            obj.print_bin("Deinterleave", codewords);
        end

        function bytes_w = dewhiten(obj, bytes)
            % dewhiten  Data Dewhitening
            %
            % input:
            %     bytes: Bytes after deinterleaving
            % output:
            %     bytes_w: Bytes after dewhitening

            len = length(bytes);
            bytes_w = bitxor(uint8(bytes(1:len)), obj.whitening_seq(1:len));
            obj.print_bin("Dewhiten", bytes_w);
        end

        % 汉明码解码
        function nibbles = hamming_decode(obj, codewords, rdd)
            % hamming_decode  Hamming Decoding
            %
            % input:
            %     codewords: Codewords after deinterleaving
            %     rdd: Bits with redundancy
            % output:
            %     nibbles: Nibbles after hamming decoding

            p1 = LoRaPHY.bit_reduce(@bitxor, codewords, [8 4 3 1]);
            p2 = LoRaPHY.bit_reduce(@bitxor, codewords, [7 4 2 1]);
            p3 = LoRaPHY.bit_reduce(@bitxor, codewords, [5 3 2 1]);
            p4 = LoRaPHY.bit_reduce(@bitxor, codewords, [5 4 3 2 1]);
            p5 = LoRaPHY.bit_reduce(@bitxor, codewords, [6 4 3 2]);
            function pf = parity_fix(p)
                switch p
                    case 3 % 011 wrong b3
                        pf = 4;
                    case 5 % 101 wrong b4
                        pf = 8;
                    case 6 % 110 wrong b1
                        pf = 1;
                    case 7 % 111 wrong b2
                        pf = 2;
                    otherwise
                        pf = 0;
                end
            end
            if obj.hamming_decoding_en
                switch rdd
                    % TODO report parity error
                    case {5, 6}
                        nibbles = mod(codewords, 16);
                    case {7, 8}
                        parity = p2*4+p3*2+p5;
                        pf = arrayfun(@parity_fix, parity);
                        codewords = bitxor(codewords, uint16(pf));
                        nibbles = mod(codewords, 16);
                    otherwise
                        % THIS CASE SHOULD NOT HAPPEN
                        error('Invalid Code Rate!');
                end
            else
                nibbles = mod(codewords, 16);
            end
            obj.print_bin("Hamming Decode", codewords);
        end

        function print_bin(obj, flag, vec, size)
            if obj.is_debug
                if nargin == 3
                    size = 8;
                end
                len = length(vec);
                fprintf("%s:\n", flag);
                for i = 1:len
                    fprintf("%s\n", dec2bin(round(vec(i)), size));
                end
                fprintf("\n");
            end
        end

        function print_hex(obj, flag, vec)
            if obj.is_debug
                len = length(vec);
                fprintf("%s: ", flag);
                for i = 1:len
                    fprintf("%s ", dec2hex(round(vec(i))));
                end
                fprintf("\n");
            end
        end

        function log(obj, flag, data)
            if obj.is_debug
                fprintf("%s: ", flag);
                len = length(data);
                for i = 1:len
                    fprintf("%d ", data(i));
                end
                fprintf("\n");
            end
        end
    end
    methods(Static)
        function b = bit_reduce(fn, w, pos)
            b = bitget(w, pos(1));
            for i = 2:length(pos)
                b = fn(b, bitget(w, pos(i)));
            end
        end

        function w = word_reduce(fn, ws)
            w = ws(1);
            for i = 2:length(ws)
                w = fn(w, ws(i));
            end
        end
    end
end

