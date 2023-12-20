%%
fre_dif_1 = zeros(1, 201);
fre_dif_2 = zeros(1, 201);
fre_dif_num = zeros(1,201);
for i = 1:size(result_array, 2)
    index = result_array(1, i) + 101;
%     disp(index);
    fre_dif_num(index) = fre_dif_num(index) + 1;
    fre_dif_1(index) = fre_dif_1(index) + result_array(2, i);
    fre_dif_2(index) = fre_dif_2(index) + result_array(3, i);
end
% index_tmp = find(fre_dif_num < 50);
% fre_dif_num(index_tmp) = 0;
subplot(2,1,1);
plot((-100:100)*500000/100, fre_dif_1./fre_dif_num);
subplot(2,1,2);
plot((-100:100)*500000/100, fre_dif_2./fre_dif_num);
    
%%
bin_num = zeros(1, 201);
amp_result = zeros(1, 201);
for i = 1:size(bin_peak, 2)
    index = bin_peak(1, i);
    bin_num(index) = bin_num(index) + 1;
    amp_result(index) = amp_result(index) + bin_peak(2, i);
end
plot(amp_result./bin_num,'.');
% plot(amp_result./bin_num);
%%
tmp = 1024;
bin_num = zeros(1, tmp*2+1);
fre_dif_1 = zeros(1, tmp*2+1);
fre_dif_2 = zeros(1, tmp*2+1);
for i = 1:size(result_array, 2)
    index = result_array(1, i) + tmp + 1;
    bin_num(index) = bin_num(index) + 1;
    fre_dif_1(index) = fre_dif_1(index) + result_array(2, i);
    fre_dif_2(index) = fre_dif_2(index) + result_array(3, i);
end
index_tmp = find(bin_num < 50);
bin_num(index_tmp) = 0;
bin_num(1:16) = 0;
bin_num(end-15:end) = 0;
bin_num(tmp-16:tmp+16) = 0;
subplot(2,1,1);
plot((-tmp:tmp)*500000/tmp, fre_dif_1./bin_num, '.');
subplot(2,1,2);
plot((-tmp:tmp)*500000/tmp, fre_dif_2./bin_num, '.');
% plot(amp_result./bin_num);
%% 统计相位
tmp = 1024;
bin_num = zeros(1, tmp*2+1);
phase_num = zeros(1, tmp*2+1);
phase_num_tmp = zeros(1, tmp*2+1);
for i = 1:size(result_array, 2)
    if result_array(1, i) < 10000
        index = result_array(1, i) + tmp + 1;
        bin_num(index) = bin_num(index) + 1;
        phase_num(index) = phase_num(index) + result_array(2, i);
%         phase_num_tmp(index) = phase_num_tmp(index) + result_array_tmp(2, i);
    end
end
index_tmp = find(bin_num < 50);
bin_num(index_tmp) = 0;
% bin_num(1:16) = 0;
% bin_num(end-15:end) = 0;
% bin_num(tmp-16:tmp+16) = 0;
plot((-tmp:tmp)*500000/tmp, phase_num./bin_num, '.b'); hold on;
% plot((-tmp:tmp)*500000/tmp, phase_num_tmp./bin_num, '.r');
%% 统计频率泄露FFT measurement
map_value = containers.Map();
map_count = containers.Map();
for i = 1:size(result_array, 2)
    tmp = num2str(result_array(1, i));
    if map_value.isKey(tmp)
        map_value(tmp) = map_value(tmp) + result_array(2, i);
        map_count(tmp) = map_count(tmp) + 1;
    else
        map_value(tmp) = result_array(2, i);
        map_count(tmp) = 1;
    end
end
tmp_result_array = zeros(2, 10); % 50e3-300, 300-475, 475-650, 650-825, 825-1000
for key = keys(map_value)
    key_tmp = str2double(key);
    if key_tmp > 0
        if key_tmp <= 300e3
            tmp_result_array(1, 6) = tmp_result_array(1, 6) + map_value(cell2mat(key));
            tmp_result_array(2, 6) = tmp_result_array(2, 6) + map_count(cell2mat(key));
        else
            index = ceil((key_tmp - 300e3)/175e3) + 6;
            tmp_result_array(1, index) = tmp_result_array(1, index) + map_value(cell2mat(key));
            tmp_result_array(2, index) = tmp_result_array(2, index) + map_count(cell2mat(key));
        end
    else
        if key_tmp >= -300e3
            tmp_result_array(1, 5) = tmp_result_array(1, 5) + map_value(cell2mat(key));
            tmp_result_array(2, 5) = tmp_result_array(2, 5) + map_count(cell2mat(key));
        else
            index = 5 - ceil((-key_tmp - 300e3)/175e3);
            tmp_result_array(1, index) = tmp_result_array(1, index) + map_value(cell2mat(key));
            tmp_result_array(2, index) = tmp_result_array(2, index) + map_count(cell2mat(key));
        end
    end
end
% subplot(2,1,1);
bar(tmp_result_array(1,:)./tmp_result_array(2,:));
set(gca,'XTickLabel', {'-1000--825','-825--650','-650--475','-475--300','-300--50','50-300','300-475','475-650','650-825','825-1000'});
% subplot(2,1,2);
% plot(result_array(1,:), result_array(3,:), '.');
%% 单统计每个信道的跳频结果
tmp_result = zeros(1,10);
tmp_result(6:10) = result_array(1, 1:2:9)./result_array(2, 1:2:9);
tmp_result(1:5) = result_array(1, 10:-2:2)./result_array(2, 10:-2:2);
bar(tmp_result);
set(gca,'XTickLabel', {'-434.175','-433.975','-433.775','-433.575','-433.375','433.375','433.575','433.775','433.975','434.175'});
            