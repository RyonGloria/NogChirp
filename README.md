# 相同带宽的非正交重叠信道解码

## Methods

1. 自下而上扫描窗口内的信号，目标 chirp 峰值应该是连续的线性增长。

2. 首先对整个目标窗口做 de-chirp，得到候选峰，滑动窗口自左向右，做差值剔除掉各个非正交的候选峰。

3. 对整个窗口做 de-chirp，得到候选峰，滤掉其它信道，再做 de-chirp，对候选峰能量进筛选。

4. 对重叠部分的窗口和非重叠部分的窗口滤波，De-chirp 后取 Bin 值的交集。

5. 滑动窗口中不同 Bin 值对应的能量方差。

## LoRa 物理层编解码设计

- 编码过程的数据形式
  1. **字符串**：'helloworldloraexp'
  2. `ASCLL`**字节(D)**：

      104  101  108  108  111  119  111  114  108  100  108  111  114   97  101  120  112

  3. `CRC` 后的**字节(D)**：

      104  101  108  108  111  119  111  114  108  100  108  111  114   97  101  120  112  0  232

  4. `白化` 后的**字节(D)**：数据白化，可以提高数据的随机性，使得数据更接近于白噪声。

      151  155  144  148  159  150  173  247  103  115   67  49   206  25  148  155  182  0  232

  5. `汉明编码`后的**字节(D)**：

      23  57  75  57  0  57  180  57  15  9  6  9  29  10  23  15  23  6  3  23  3  20  17  3  30  12  9  17  20  9  27  9  6  27  0  0  24  30

  6. `交织编码`后的**符号(D)**：

      175  130  80  213  190  87  16  128  475  1010  501  407  21  823  563  562  56  677  46  282  388  869  370

  7. 交织编码后`格雷码`对应的**符号(D)**：

      809  1009  385  613  849  405  125  1021  366  676  346  283  26  550  990  989  48  826  53  493  264  583  420

- 编码率、数据长度、扩频因子之间的推算

    > *plen*: payload length
    >
    > cr: code rate (1: $\frac45$; 2: $\frac46$; 3: $\frac47$; 4: $\frac48$​;)
    >
    > sf: spread factor

$$
\frac{4}{4 + cr} = \frac{plen \times 8}{2 \times plen \times 8} \times           \frac{sf \times 8 } {(4 + cr) \times sf}
$$
$$
Symbols = \frac{4 + cr}{ 2 \times plen \times sf }
$$



## 实验设计

> **Channel 1**：
>
> - 中心频率：433.46875 MHz
> 
>**Channel 2**：
> 
>- 中心频率：433.5 MHz
> 
> **Channel 3**：
>
> - 中心频率：433.53125 MHz

![NogChannel](Figure/NogChannel.png)

### Payload

1. - 发送端：

       ```c++
       String data = "helloworldloraexp";
       radio.implicitHeader(17)
       radio.setCRC(false, false)
       // bytes: 104 101 108 108 111 119 111 114 108 100 108 111 114 97 101 120 112
       // bins: 809 1009 385 613 849 405 125 1021 366 676 346 283 26 550 990 989 48 826 821 84 264 586 420
       ```
   
2. - 发送端：

       ```c++
       String data = "helloworldloraexp";
       radio.implicitHeader(17)
       radio.setCRC(true, false)
       // bytes: 104 101 108 108 111 119 111 114 108 100 108 111 114 97 101 120 112
       // bins: 809 1009 385 613 849 405 125 1021 366 676 346 283 26 550 990 989 48 826 53 493 264 583 420
       ```
   
   - 接收端：
   
     ```matlab
     phy = LoRaPhyLay(sf, bw, samplesRate);
     phy.has_header = 0;         % implicit header mode
     phy.crc = 1;                % enable payload CRC checksum
     phy.payload_len = 17;       % payload length 17 bytes
     ```
   
     

### 文件目录

- 虚拟机共享文件夹：\\192.168.3.102\e\share\samples\

### Payload Detection

1. CIC 方法：先检测 SFD，根据 SFD 位置检测前面的 preamble（低信噪比下检测是否会出问题）
2. 

### 实验记录

CIC 解码结果：

 CH1 payload1: 失败（该信道能量较低，可能被直接滤掉了）

 CH2 payload1: [810 722 386 614 850 593 126 1022 842 677 347 862 581 551 574 113 806 827 822 85 421 587 421]

 CH3 payload2: [810 1010 386 614 850 406 126 1022 367 677 347 284 27 551 1004 1006 49 836 822 85 265 587 421]

准确率：0/23 14/23 20/23  

1. - **问题**：提出一个新的解决方案后，直接通过实验验证难度大，可行性低
   - **答**：编写单 Chirp 信号仿真平台，包括 downchirp、baseupchirp、冲突 chirp 和非正交冲突 chirp 的生成，在此平台上进行方法的可行性验证。

2. - **问题**：实验中发现不同硬件解码的 Bin 值不一致
   - **答**：nf95 由于接口配置问题，在实际使用中不同硬件解码的 Bin 值不一致，改换 Radiolib 库并且在每次 setup 阶段 reset 对应的接口。

3. - **问题**：难以控制不同节点发送数据产生冲突
   - **答**：编写发送端程序，通过 1 个节点控制 3 个节点分别在 Channel 1、2、3 发送数据。

4. - **问题**：实现一控多之后，发现节点几乎同时发送数据
   - **答**：Channel 1 发送数据时间偏移 50 ms，Channel 3发送数据时间偏移 100 ms

5. - **问题**：量哥指出直接进行 3 个非正交信道的解码工作实现难度大
   - **答**：可以先考虑两个非正交信道的解码，例如原信道 Channel 1 占带宽 125 KHz，此时增加 50 KHz 带宽，即可划分两个非正交信道 Channel 1 和 2。即通过增加少量的带宽资源可实现更多信道的划分。

6. - **问题**：互相关性能否用在 SFD 的检测中？
   - **答**：参考 CIC 方法

7. - **问题**：由于信号是含正负方向的，滤波器无法滤掉指定的频率
   - **答**：将信号乘上  $e^{j 2\pi t f_s}$ 其中 $f_s$ 是偏移的频率，之后再进行滤波。

8. - **问题**：连续线性变化计算开销大，是否可以从频率角度直接考虑两个窗口的交集
   - **答**：对重叠部分的窗口和非重叠部分的窗口滤波，De-chirp 后取 Bin 值的交集。

9. - **问题**：切分的滑动窗口，de-chirp 结果存在峰值偏移的问题
   - **答**：是否可以通过设置置信区间解决该问题

10. - **问题**：（240112 量哥）可以用三种 downchirp：1.只在信道 1 未重叠部分的 downchirp；2. 只在信道2未重叠的 downchirp； 3. 只在信道 1 和信道 2 重叠部分的 downchirp
    - **答**： 不可行；过滤掉频谱的

11. 重要：*网关侧采集不同设备发送的 payload 会产生不同的 CFO/FFO，非正交重叠信道本质上就是频率偏移。*
    - CH1-2-3 CIC 解：能解出 CH2-3，CH1 因为能量较低被过滤
    - CH2-3-1 CIC 解：能解出 CH2-3，CH1 检测不到

12. 增加少些信道资源，能指数级增加吞吐量

13. 非周期性截断导致频谱泄露，可通过加窗的方式缓解

14. 即使分离不同信道信号，同信道信号解码仍有很大难度

15. - **问题**：由于旁瓣效应导致峰值偏移只会是 Bin 值偏移 1~2，可以考虑只对最后两位进行编码。

    - **答**： 

## 函数说明

### LoRaPhyLay.m

> 该函数作用是将 Symbols 对应的 Bin 值转为 Bytes，具体实例如下

```matlab
sf = 10;                    % spreading factor SF10
bw = 125e3;                 % bandwidth 125 kHz
sample_rate = 2e6;          % sampling rate 2 MHz

phy = LoRaPhyLay(sf, bw, sample_rate);
phy.has_header = 0;         % implicit header mode
phy.cr = 1;                 % code rate = 4/8 (1:4/5 2:4/6 3:4/7 4:4/8)
phy.crc = 1;                % enable payload CRC checksum
phy.payload_len = 17;       % payload length 17 bytes
phy.is_debug = 0;           % enable debug mode

str = 'helloworldloraexp';
ascii_values = double(str);
fprintf("[encode] data:\n");
disp(ascii_values);

% Encoding
symbols = phy.encode(ascii_values');
fprintf("[encode] symbols:\n");
disp(symbols');

% Decoding
[data, checksum] = phy.decode(symbols);
fprintf("[decode] data:\n");
disp(data');

fprintf("[decode] checksum:\n");
disp(checksum');

fprintf("[decode] data without checksum:\n");
disp(data(1:17)');

data_str = char(data(1:17)');
disp(data_str);
```

