> 该文档描述了目录下所有函数及相关代码的释义
>
> $BW$：带宽
> $SF$：扩频因子
> $CR$：编码率
> $T_s$：Symbol 持续时间
> $N_s = \frac{1}{T_s}$：表示为 1s 内 Symbol 个数

首先，我们可以写出 Chirp 信号瞬时频率与时间 $t$ 的关系，即
$$
f_{chirp}(t)=-\frac{BW}{2}+\frac{BW}{T_s}t
$$
实际上，对于一个复函数 $y(t) = e^{j\phi(t)}$ 为虚数单位，即 -1 的平方根，其瞬时频率的表达式为
$$
f_y(t) = \frac{\phi'(t)}{2\pi}
$$
因此了得到 Chirp 信号在时域上的表达式，我们可以求出 $\phi(t)$，即

$$
\begin{aligned}
\phi(t)
& = \int f_{chirp}(t) dt \\
& = 2\pi t (-\frac{BW}{2} + \frac{BW}{2T_s} t)
\end{aligned}
$$
于是我们就得到了 Chirp 在时域上复信号
$$
y_{upchirp}(t) = e^{j2\pi t (-\frac{BW}{2} + \frac{BW}{2T_s} t)}
$$
我们称上述 Chirp 信号为一个标准 upchirp，类似的我们可以构造标准 downchirp，写作

$$
y_{downchirp}(t) = e^{j2\pi t (\frac{BW}{2} - \frac{BW}{2T_s} t)}
$$
如果我们将标准 upchirp 与标准 downchirp 做乘法，即
$$
\begin{aligned}
s(t) & = y_{upchirp}(t)y_{downchirp}(t) \\
& = e^{j2\pi t (-\frac{BW}{2} + \frac{BW}{2T_s} t)}e^{j2\pi t (\frac{BW}{2} - \frac{BW}{2T_s} t)} \\ 
& = e^0 \\
& = 1
\end{aligned}
$$
有时我们将 upchirp 写作

$$
C = e^{j2\pi t (-\frac{BW}{2} + \frac{BW}{2T_s} t)}
$$
而 downchirp 写作

$$
C^{-1} = e^{j2\pi t (\frac{BW}{2} - \frac{BW}{2T} t)}
$$
这符合大多数公理系统的逆运算，即满足 $CC^{-1} = 1$。

## 函数 buildIdealchirp

### 代码：

```matlab
function obj = buildIdealchirp(obj, f_temp)
            cmx = 1+1*1i;
            pre_dir = 2*pi;
            f0 = f_temp + obj.loraSet.bw/2;                           % 设置理想upchirp和downchirp的初始频率
            f1 = -f_temp+obj.loraSet.bw/2;
            d_symbols_per_second = obj.loraSet.bw / obj.loraSet.fft_x;
            T = -0.5 * obj.loraSet.bw * d_symbols_per_second;
            d_samples_per_second = obj.loraSet.sample_rate;        % sdr-rtl的采样率
            d_dt = 1/d_samples_per_second;         % 采样点间间隔的时间
            t = d_dt*(0:1:obj.loraSet.dine-1);
            % 计算理想downchirp和upchirp存入d_downchirp和d_upchirp数组中（复数形式）
            obj.idealDownchirp = cmx * (cos(pre_dir .* t .* (f0 + T * t)) + sin(pre_dir .* t .* (f0 + T * t))*1i);
            obj.idealUpchirp = cmx * (cos(pre_dir .* t .* (f1 + T * t) * -1) + sin(pre_dir .* t .* (f1 + T * t) * -1)*1i);
        end
```

### 释义:

```matlab
f0 = f_temp + obj.loraSet.bw/2;
f1 = -f_temp+obj.loraSet.bw/2;
```


$$
f_0 = f_1 = \frac{BW}{2}
$$

其中，$f_{temp}$ 指频率（BW）上下的偏移

```matlab
d_symbols_per_second = obj.loraSet.bw / obj.loraSet.fft_x;
```

fft_x 为 $2^{SF}$，此处公式表示为
$$
N_s = \frac{BW}{2^{SF}}
$$
```matlab
T = -0.5 * obj.loraSet.bw * d_symbols_per_second;
```

公式表示为
$$
T = \frac{-BW*N_s}{2} = \frac{-BW}{2T_s}
$$

```matlab
obj.idealDownchirp = cmx * (cos(pre_dir .* t .* (f0 + T * t)) + sin(pre_dir .* t .* (f0 + T * t))*1i);
```

公式表示为
$$
\begin{aligned}
y_{downchirp}(t) & = (1+i)\bigg(\cos{\big(2\pi t(f_0+T \times t)\big)}+\sin{\big(2\pi t(f_0+T \times t)\big)\times i\bigg)} \\
& = \sqrt{2} \cos(x+\frac{\pi}{4}) + i \times \sqrt{2} \sin(x+\frac{\pi}{4})  , x = 2\pi t(\frac{BW}{2} - \frac{BW}{2T_s} t)\\
& = \sqrt{2} e^{i2\pi t (\frac{BW}{2} - \frac{BW}{2T_s} t + \frac{\pi}{4})}
\end{aligned}
$$


## 
