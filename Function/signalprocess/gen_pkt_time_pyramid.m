function [T]=gen_pkt_time_pyramid(node_num,pkt_rate,pkt_num,pkt_lengthT,fs)
    % 生成泊松分布的发包时间
    % node_num: 节点数目
    % pkt_rate: 发包速率 每分钟发包数目
    % pkt_num: 发包数目 一共发包次数
    % pkt_length: 包长，控制发包间隔大于一个数据包长度
    % fs: 采样率
    lambda = pkt_rate/(60*fs);%节点数据包产生速率，泊松到达率
    u = unifrnd(0,1,node_num,pkt_num);
    t = -log(u)/lambda+pkt_lengthT; %发包时间间隔符合指数分布
    t = [zeros(node_num,1),t];
    T = zeros(node_num,pkt_num);%发包时间
    for i = 2 : pkt_num
        T(:,i) = sum(t(:,1:i),2);
    end