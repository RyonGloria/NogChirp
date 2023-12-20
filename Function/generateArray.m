function output_matrix = generateArray(num, range)
    % 参数
    mask = uint64(2^31);  
    multiplier = uint64(1103515245); 
    c = uint64(12345);  
    
    % 确保输入在0到1023之间
    x = uint64(mod(num, 1024));
    
    % 生成1x1000的矩阵，元素在0到15之间
    output_matrix = uint64(zeros(1, 1000));
    for i = 1:1000
        x = mod(multiplier * x + c, mask);
        % 将 uint64 转换为字符串
        data_str = num2str(x);
        data_str = fliplr(data_str);
        
        % 获取第五位数的字符
        if length(data_str) < 5
            fifth_digit = 0;
        else
            fifth_digit = data_str(5);
        end
        
        % 转换字符为数字
        fifth_digit_number = str2double(fifth_digit);

        output_matrix(i) = uint64(mod(fifth_digit_number, range) + 1);
    end
end


% array = [];
% for i = 0:1023
%     array = [array; generateArray(i, 4)];
% end
