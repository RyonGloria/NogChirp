function accuracy = calculateAccuracy(array1, array2)
    % 获取两个数组的长度
    len1 = length(array1);
    len2 = length(array2);

    % 如果array2的长度小于array1，则将其补齐为array1的长度
    if len2 < len1
        array2 = [array2, zeros(1, len1 - len2)];
    elseif len2 > len1
        array2 = array2(1:len1);
    end

    % 计算正确匹配的元素数量
    correctMatches = sum(array1 == array2);

    % 计算正确率
    accuracy = correctMatches / len1;
end
