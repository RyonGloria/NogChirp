function [resultSort] = sortFileByTime(FileArray)
    fileTime = zeros(size(FileArray, 2), size(FileArray, 1));
    for i = 1:size(FileArray, 1)
        for row = 1:size(FileArray, 2)
            fileTime(row, i) = FileArray(i, row).datenum;
        end
    end
    resultSort = zeros(size(FileArray, 2), size(FileArray, 1));
    for i = 1:size(FileArray, 2)
        [~, resultSort(i, :)] = sort(fileTime(i, :));
    end