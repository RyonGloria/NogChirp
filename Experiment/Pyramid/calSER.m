fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器
bin_ref = load('E:\Pyramid_samples\bin_ref.txt')';
bin = textread('E:\Pyramid_samples\bin.txt','','delimiter',',','emptyvalue',NaN);
if size(bin,2) < length(bin_ref)
    bin = [bin(), zeros(size(bin), length(bin_ref) - size(bin,2))];
end

pkg_num = 40;
result = 0;
for i = 1:size(bin, 1)
    bin_tmp = bin(i, 1:length(bin_ref));
    tmp = (bin_tmp - (bin_ref-1));
    result = result + sum((bin_tmp - (bin_ref-1)) == 0);
end
num_all = pkg_num * length(bin_ref);
fprintf('SRR is %f\n', result/num_all);
toc;
fclose all;