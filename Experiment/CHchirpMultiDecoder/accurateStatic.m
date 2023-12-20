%% 统计binTrueRate
% collNum = 8;
% binNum = 20;
% % resultArr = zeros(1, collNum*size(record, 2));
% resultArr = zeros(1, binNum+1);
% for index = 1:size(record, 2)
%     binRecord = record{12, index};
%     for pkgIndex = 1:collNum
%         if pkgIndex < length(binRecord)
%             trueBinNum = binRecord(pkgIndex) * binNum;
%             resultArr(trueBinNum+1) = resultArr(trueBinNum+1) + 1;
%         else
%             resultArr(1) = resultArr(1) + 1;
%         end
%     end
% end
% resultArr = resultArr/sum(resultArr);
% bar(resultArr);

%% 统计binOffset
collNum = 8;
fft_x = 1024;
sf = 10;
true_bin = importdata(strcat('.\Code\Config\bin\SF', string(sf), '.txt'))';
% resultArr = zeros(1, collNum*size(record, 2));
resultArr = zeros(1, fft_x*2 + 1);
for index = 1:size(record, 2)
    binRecord = record{10, index};
    for pkgIndex = 1:collNum
        if pkgIndex < length(binRecord)
            binArr = binRecord{pkgIndex};
            for binIndex = 1:length(true_bin)
                resultArr(binArr(binIndex)-true_bin(binIndex)+fft_x+1) = resultArr(binArr(binIndex)-true_bin(binIndex)+fft_x+1) + 1;
            end
        end
    end
end
% resultArr = resultArr/sum(resultArr);
resultArr(1025) = 0;
bar(resultArr);