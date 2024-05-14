function errorDist = calculateBitError(array1, array2, sf)
    % 获取两个数组的长度
    len1 = length(array1);
    len2 = length(array2);

    % 如果array2的长度小于array1，则将其补齐为array1的长度
    if len2 < len1
        array2 = [array2, zeros(1, len1 - len2)];
    elseif len2 > len1
        array2 = array2(1:len1);
    end

    % errorIndices = find(array1 ~= array2);
    % Find indices where elements differ, excluding NaN values
    validIndices = find(~isnan(array1) & ~isnan(array2));
    errorIndices = validIndices(array1(validIndices) ~= array2(validIndices));
    if ~isempty(errorIndices)
        DecTrue = array1(errorIndices);
        DecErr = array2(errorIndices);

        binDist = bitxor(DecTrue, DecErr);
        valDist = min(abs(abs(DecErr - DecTrue) - 2^sf), abs(DecTrue - DecErr));

        errorDist = [valDist; binDist];
        % errorPositions = find(dec2bin(errorBits, 10) == '1');
        disp(['Errors at indices: ', num2str(errorIndices), ' Orginal: ', num2str(DecTrue), ' Wrong: ', num2str(DecErr)]);
    else
        errorDist = [];
    end
end
