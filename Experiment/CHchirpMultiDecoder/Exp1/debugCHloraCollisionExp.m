%% 统计解码出包的数目
recordLength = size(record, 2);
colNumArr = zeros(1, recordLength);
decodePkgNumArr = zeros(1, recordLength);
for index = 1:recordLength
    colNumArr(index) = record{1, index};
    decodePkgNumArr(index) = length(record{10, index});
end

collEqual = zeros(1, record{1, end});
collSmall = zeros(1, record{1, end});
collBig = zeros(1, record{1, end});
for index = 1:recordLength
    if decodePkgNumArr(index) == colNumArr(index)
        collEqual(colNumArr(index)) = collEqual(colNumArr(index)) + 1;
    elseif decodePkgNumArr(index) > colNumArr(index)
        collBig(colNumArr(index)) = collBig(colNumArr(index)) + 1;
    else
        collSmall(colNumArr(index)) = collSmall(colNumArr(index)) + 1;
    end
end