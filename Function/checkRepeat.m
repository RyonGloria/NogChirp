% 生成一个随机的1000x1024的矩阵（假设为matrix）
% matrix = randi([0, 1], 1000, 1024);  % 随机生成一个二值矩阵

% 使用unique函数查看是否有重复行
[unique_matrix, ~, idx] = unique(array, 'rows', 'stable');

% 如果有重复行，unique_matrix 将包含唯一的行，idx 将包含索引，以指示每行在unique_matrix中的位置

if size(unique_matrix, 1) < size(array, 1)
    fprintf('存在完全相同的两行，它们位于行 %d 和行 %d。\n', idx(1), idx(2));
else
    fprintf('矩阵中没有完全相同的两行。\n');
end
