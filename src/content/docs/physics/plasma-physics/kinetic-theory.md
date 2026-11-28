---
title: 運動論
description: 微視的な位相空間の保存則から Vlasov 方程式を導き、共鳴と初期値応答を解析する。
notation:
  - symbol: 'f_s'
    meaning: 一粒子分布
  - symbol: '\zeta'
    meaning: 位相速度を熱速度で規格化した量
  - symbol: 'Z'
    meaning: プラズマ分散関数
  - symbol: '\varepsilon'
    meaning: 縦誘電関数
---

## この章の目標

- 粒子分布の平均と相関を区別する。
- 初期値問題を Laplace 変換で解く。
- Landau 減衰と捕捉粒子の定常解を説明する。

## 分布関数と Klimontovich 方程式

一つの実現の粒子位置と速度を全て記録するには、位相空間に点測度を置く。種 $s$ の粒子数 $N_s$ に対して
$$
\mathcal F_s(\vb x,\vb v,t)=\sum_{a=1}^{N_s}\delta^{(3)}(\vb x-\vb x_a(t))\delta^{(3)}(\vb v-\vb v_a(t))
$$
とする。これを :term[Klimontovich 分布]{en="Klimontovich distribution" yomi="くりもんとびっちぶんぷ"} と呼ぶ。粒子の自己力は除き、点力での積は各軌道上の力として解釈する。滑らかな試験関数 $g$ を掛けて積分すると $\sum_a g(\vb x_a,\vb v_a)$ である。時間微分は $\sum_a(\vb v_a\cdot\nabla_x g+\vb a_a\cdot\nabla_v g)$ になるので、分布の意味で
$$
\partial_t\mathcal F_s+\nabla_x\cdot(\vb v\mathcal F_s)+\nabla_v\cdot(\vb a_s\mathcal F_s)=0
\label{eq-plasma-klimontovich}
$$
である。Lorentz 加速度 $\vb a_s=q_s(\vb E+\vb v\times\vb B)/m_s$ は速度発散がゼロであり、流れは位相体積を保存する。

:::example[自由飛行]
場がないと $\vb x_a(t)=\vb x_a(0)+t\vb v_a(0)$ である。よって $\mathcal F_s(\vb x,\vb v,t)=\mathcal F_s(\vb x-t\vb v,\vb v,0)$ となる。
:::

:::example[粗視化の違い]
同一の速度 $v_0$ の粒子二個が $x_1,x_2$ にいる一次元模型では、密度の Fourier 成分は $e^{-ikx_1}+e^{-ikx_2}$ である。$k(x_2-x_1)=\pi$ ならその成分はゼロである。全粒子数はいずれの位相差でも二個であり、密度の積分と個々の情報を区別できる。
:::

### 自由飛行の粗視化で密度が減衰する理由

一次元の初期分布を $f(x,v,0)=F_M(v)[1+a\cos(kx)]$ とする。$F_M$ は平均 $U$、分散 $v_t^2$ の Maxwell 分布である。力がなければ特性をたどり $x$ を $x-vt$ に置き換える。従って空間波の振幅は速度の Fourier 変換 $\int F_M(v)e^{-ikvt}\dd v=n_0e^{-ikUt-(kv_tt)^2/2}$ で決まる。Gaussian の変換は、波数について微分して速度について部分積分すれば $\hat F_M'=(-iU-v_t^2\tau)\hat F_M$ となることから得られる。密度振幅は消えても速度ごとの相対振幅は $a$ のままである。この例は自己電場を外した位相混合であり、自己無撞着な Landau 極を計算したものではない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $a$ | 初期密度波の振幅 | 1 |
| $s$ | 規格化時間 $kv_t t$ | 1 |
| $u$ | 平均流と熱速さの比 $U/v_t$ | 1 |
| $\theta$ | 観測点の空間位相 $kx$ | 1 |

速度積分で残る包絡を $D$ とする。Gaussian の Fourier 変換である。

$$
D=e^{- \frac{s^{2}}{2}}
$$

観測点の合成位相を $\Phi$ とする。平均流は波を移流させる。

$$
\Phi=\theta - s u
$$

局所の密度比を $n/n_0$ とする。速度分布の面積を取る。

$$
n/n_0=D a \cos{\left(\Phi \right)} + 1
$$

密度の Fourier 振幅を $A_n$ とする。位置を変えて最大振幅を読む。

$$
A_n=D a
$$

速度方向に一周する間隔を $\Delta v/v_t$ とする。時間とともに細かい構造になる。

$$
\Delta v/v_t=\frac{2 \pi}{s}
$$

初期密度振幅に対する比を $R$ とする。速度分布自体の粗視化とは区別する。

$$
R=\frac{A_{n}}{a}
$$

:::example[初期の緩やかな減衰]
振幅と速度空間の間隔を求める。

入力値を次のように置く。

$$
\begin{aligned}
a&=\frac{1}{10} \\s&=\frac{1}{2} \\u&=0 \\\theta&=0
\end{aligned}
$$

これらを順に代入する。

速度積分で残る包絡は

$$
D=e^{- \frac{\left(\frac{1}{2}\right)^{2}}{2}}\simeq 0.8825
$$

である。

観測点の合成位相は

$$
\Phi=- \left(0\right) \left(\frac{1}{2}\right) + \left(0\right)=0
$$

である。

局所の密度比は

$$
n/n_0=\left(\frac{1}{10}\right) \left(e^{- \frac{1}{8}}\right) \cos{\left(\left(0\right) \right)} + 1\simeq 1.0882
$$

である。

密度の Fourier 振幅は

$$
A_n=\left(\frac{1}{10}\right) \left(e^{- \frac{1}{8}}\right)\simeq 0.08825
$$

である。

速度方向に一周する間隔は

$$
\Delta v/v_t=\frac{2 \pi}{\left(\frac{1}{2}\right)}\simeq 12.566
$$

である。

初期密度振幅に対する比は

$$
R=\frac{\left(\frac{1}{10 e^{\frac{1}{8}}}\right)}{\left(\frac{1}{10}\right)}\simeq 0.8825
$$

である。

幅の大きな速度構造があり、積分後の減衰は弱い。
:::

:::example[熱速さによる位相混合]
平均流ゼロでも密度が減るか確かめる。

入力値を次のように置く。

$$
\begin{aligned}
a&=\frac{1}{10} \\s&=2 \\u&=0 \\\theta&=0
\end{aligned}
$$

これらを順に代入する。

速度積分で残る包絡は

$$
D=e^{- \frac{\left(2\right)^{2}}{2}}\simeq 0.13534
$$

である。

観測点の合成位相は

$$
\Phi=- \left(0\right) \left(2\right) + \left(0\right)=0
$$

である。

局所の密度比は

$$
n/n_0=\left(\frac{1}{10}\right) \left(e^{-2}\right) \cos{\left(\left(0\right) \right)} + 1\simeq 1.0135
$$

である。

密度の Fourier 振幅は

$$
A_n=\left(\frac{1}{10}\right) \left(e^{-2}\right)\simeq 0.013534
$$

である。

速度方向に一周する間隔は

$$
\Delta v/v_t=\frac{2 \pi}{\left(2\right)}\simeq 3.1416
$$

である。

初期密度振幅に対する比は

$$
R=\frac{\left(\frac{1}{10 e^{2}}\right)}{\left(\frac{1}{10}\right)}\simeq 0.13534
$$

である。

速度階級ごとの空間位相がずれ、密度で相殺する。
:::

:::example[並進を伴う速度分布]
平均流を追う観測点で密度を求める。

入力値を次のように置く。

$$
\begin{aligned}
a&=\frac{1}{10} \\s&=1 \\u&=2 \\\theta&=2
\end{aligned}
$$

これらを順に代入する。

速度積分で残る包絡は

$$
D=e^{- \frac{\left(1\right)^{2}}{2}}\simeq 0.60653
$$

である。

観測点の合成位相は

$$
\Phi=- \left(1\right) \left(2\right) + \left(2\right)=0
$$

である。

局所の密度比は

$$
n/n_0=\left(\frac{1}{10}\right) \left(e^{- \frac{1}{2}}\right) \cos{\left(\left(0\right) \right)} + 1\simeq 1.0607
$$

である。

密度の Fourier 振幅は

$$
A_n=\left(\frac{1}{10}\right) \left(e^{- \frac{1}{2}}\right)\simeq 0.060653
$$

である。

速度方向に一周する間隔は

$$
\Delta v/v_t=\frac{2 \pi}{\left(1\right)}\simeq 6.2832
$$

である。

初期密度振幅に対する比は

$$
R=\frac{\left(\frac{1}{10 e^{\frac{1}{2}}}\right)}{\left(\frac{1}{10}\right)}\simeq 0.60653
$$

である。

位相のずれを追跡しても熱分散から生じる減衰は残る。
:::

:::example[密度の節での観測]
ある位置の密度変化から全振幅を推定できるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=\frac{1}{10} \\s&=1 \\u&=1 \\\theta&=1 + \frac{\pi}{2}
\end{aligned}
$$

これらを順に代入する。

速度積分で残る包絡は

$$
D=e^{- \frac{\left(1\right)^{2}}{2}}\simeq 0.60653
$$

である。

観測点の合成位相は

$$
\Phi=\left(1 + \frac{\pi}{2}\right) - \left(1\right)^{2}\simeq 1.5708
$$

である。

局所の密度比は

$$
n/n_0=\left(\frac{1}{10}\right) \left(e^{- \frac{1}{2}}\right) \cos{\left(\left(\frac{\pi}{2}\right) \right)} + 1=1
$$

である。

密度の Fourier 振幅は

$$
A_n=\left(\frac{1}{10}\right) \left(e^{- \frac{1}{2}}\right)\simeq 0.060653
$$

である。

速度方向に一周する間隔は

$$
\Delta v/v_t=\frac{2 \pi}{\left(1\right)}\simeq 6.2832
$$

である。

初期密度振幅に対する比は

$$
R=\frac{\left(\frac{1}{10 e^{\frac{1}{2}}}\right)}{\left(\frac{1}{10}\right)}\simeq 0.60653
$$

である。

局所の密度は背景に等しいが、別の位置には波がある。
:::


## BBGKY 階層と Vlasov 方程式

粒子系の確率密度を縮約すると、一粒子の方程式に二粒子分布が現れる。同種 $N$ 個、粒子対の加速度 $\vb a_{ij}$ を考え、全確率密度を $P_N$、$f_j=N!/(N-j)!\int P_N\dd z_{j+1}\cdots\dd z_N$、$z_i=(\vb x_i,\vb v_i)$ と規格化する。Liouville 方程式を積分すると、積分された粒子についての発散項は消える。残るのは、数えた粒子の移流と内部の相互作用、および外の一粒子との相互作用である。置換対称性から外側の $N-j$ 項は同一で、この係数は $f_{j+1}$ の規格化に吸収される。よって
$$
\left[\partial_t+\sum_{i=1}^j\vb v_i\cdot\nabla_{x_i}
+\sum_{i\ne l\le j}\vb a_{il}\cdot\nabla_{v_i}\right]f_j
=-\sum_{i=1}^j\nabla_{v_i}\cdot\int \vb a_{i,j+1} f_{j+1}\dd z_{j+1}.
$$
これが :term[BBGKY 階層]{en="BBGKY hierarchy" yomi="びーびーじーけーわいかいそう"} である。$f_2=(1-1/N)f_1f_1+g_2$ と分ける。$g_2$ の寄与が小さい時間で捨て、大きな $N$ で自己寄与の差も落とすと場が $f_1$ から決まる。多種での平均場と Maxwell 方程式による表現は
$$
\partial_t f_s+\vb v\cdot\nabla_x f_s+
\frac{q_s}{m_s}(\vb E+\vb v\times\vb B)\cdot\nabla_v f_s=0 .
\label{eq-plasma-vlasov}
$$
これは :term[Vlasov 方程式]{en="Vlasov equation" yomi="ぶらそふほうていしき"} である。場を粒子数密度と電流から決める点が試験粒子の輸送式と異なる。$N_D\gg1$ は平均場を使う目安であり、長い時間での相関の蓄積まで無視できるとは限らない。

:::example[位相空間での二乗積分]
境界流出がないとき式に $2f_s$ を掛けると、$f_s^2$ も発散形で移流する。従って $\int f_s^2\dd^3x\dd^3v$ は保存される。
:::

### 微視的な方程式の平均で残る項

Klimontovich 分布を $\mathcal F_s=f_s+\delta f_s$、微視的な Lorentz 力を $\vb F_s=\overline{\vb F}_s+\delta\vb F_s$ と分ける。ここで $f_s=\langle\mathcal F_s\rangle$、$\overline{\vb F}_s=\langle\vb F_s\rangle$ であり、平均からのずれはそれぞれ平均ゼロである。位相空間の微視的な保存式を平均すると、力の項の積には二つの寄与が残る。従って一粒子分布は
$$
\partial_t f_s+\vb v\cdot\nabla_x f_s+
\frac{\overline{\vb F}_s}{m_s}\cdot\nabla_v f_s
=-\frac{1}{m_s}\nabla_v\cdot\langle\delta\vb F_s\,\delta f_s\rangle
$$
を満たす。ここでは平均力の速度発散がゼロであることを使った。Maxwell 方程式が場について線形なので、平均場の源は $f_s$ の電荷と電流のモーメントである。右辺は揺らぐ場と分布の相関であり、平均ゼロの揺らぎ同士でも残る。静電の場合、自己力を除いた右辺は BBGKY の二体相関の寄与に対応する。平均場近似で捨てるのは、この右辺の相関項である。

### Vlasov–Poisson の保存量

周期境界、時間に依存しない一様な中和背景と、周期電位から得る静電場を考える。式 :eqref[eq-plasma-vlasov] に $m_s v^2/2$ を掛けて積分すると粒子エネルギーの変化率は $\int\vb j\cdot\vb E\dd V$ である。Poisson 方程式を時間微分して $\varepsilon_0\partial_t\nabla\cdot\vb E=-\nabla\cdot\vb j$ とし、$\phi$ を掛けて積分する。$\vb E=-\nabla\phi$ と周期境界から $\varepsilon_0\int\vb E\cdot\partial_t\vb E\dd V=-\int\vb j\cdot\vb E\dd V$ となる。従って粒子と電場の総エネルギーは保存される。滑らかな関数 $C$ について $C(f_s)$ も位相空間で移流される。速度方向で可積分かつ無限遠の境界流束が消えれば、その位相空間積分も保存される。これらの積分が保存されることと、速度を積分した低次モーメントが同じ値を保つことは別である。

:::example[二つの解像度]
自由移流の初期摂動が $\cos kx$ と滑らかな速度依存を持てば、時刻 $t$ の位相は $k(x-vt)$ である。固定の $x$ で見ると速度の尺度は $2\pi/(kt)$ に細かくなる。$t$ を四倍にすれば速度格子に必要な間隔も四分の一にする必要がある。
:::

### 平均場近似の前後で二体相関を区別する

箱の中に $N$ 個の同種粒子を取り、粒子密度の相関を確率として扱う。$p_1,p_2$ をそれぞれ積分が一の周辺確率密度とし、$p_2=p_1p_1+g_2$ と分ける。一体方程式の左辺の力項は $\vb a_{\rm mf}\cdot\nabla_{\vb v_1}p_1$ と $(N-1)\nabla_{\vb v_1}\cdot\int\vb a_{12}g_2\dd z_2$ の和になる。ただし $\vb a_{\rm mf}=(N-1)\int\vb a_{12}p_1(z_2)\dd z_2$ であり、$\vb a_{12}$ は粒子 2 による粒子 1 の加速度である。分ける理由は $\int g_2(1,2)\dd2=0$ でも力を掛けた積分は残るためである。単純な統計的比較として、ある小胞の占有確率を $p$ とし粒子が独立なら、個数 $M=\sum_{j=1}^N I_j$ で $\mathbb E M=Np,\ \operatorname{Var}M=Np(1-p)$ となる。指示変数同士に共分散 $c$ があれば分散へ $N(N-1)c$ が加わる。大きい $N$ でも弱い個別相関が和に残り得る。この計算は一つの小胞に対する統計であり、Coulomb 力の特異性を含んだ Vlasov 極限の数学的証明ではない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $N$ | 全粒子数 | 1 |
| $p$ | 一粒子が小胞を占める確率 | 1 |
| $c$ | 異なる二粒子の占有指示変数の共分散 | 1 |

平均占有数を $\mathbb E M$ とする。共分散の有無に依存しない。

$$
\mathbb E M=N p
$$

対角分散を $V_0$ とする。各指示変数の分散を足す。

$$
V_0=N p \left(1 - p\right)
$$

二体相関の分散寄与を $V_2$ とする。順序付きの相異なる対を数える。

$$
V_2=N c \left(N - 1\right)
$$

全占有数の分散を $\operatorname{Var}M$ とする。二体共分散は負にもなり得る。

$$
\operatorname{Var}M=V_{0} + V_{2}
$$

密度の相対ゆらぎを $R_{\rm rms}$ とする。少数の体積では大きくなる。

$$
R_{\rm rms}=\frac{\sqrt{\operatorname{Var}M}}{\mathbb E M}
$$

独立粒子の尺度からの補正を $V_2/V_0$ とする。値が一に近いと独立性を使った推定は変わる。

$$
V_2/V_0=\frac{V_{2}}{V_{0}}
$$

:::example[独立な粗視化]
百個程度を含む小胞で相対ゆらぎを求める。

入力値を次のように置く。

$$
\begin{aligned}
N&=10000 \\p&=\frac{1}{100} \\c&=0
\end{aligned}
$$

これらを順に代入する。

平均占有数は

$$
\mathbb E M=\left(10000\right) \left(\frac{1}{100}\right)=100
$$

である。

対角分散は

$$
V_0=\left(10000\right) \left(\frac{1}{100}\right) \left(1 - \left(\frac{1}{100}\right)\right)=99
$$

である。

二体相関の分散寄与は

$$
V_2=\left(0\right) \left(10000\right) \left(\left(10000\right) - 1\right)=0
$$

である。

全占有数の分散は

$$
\operatorname{Var}M=\left(0\right) + \left(99\right)=99
$$

である。

密度の相対ゆらぎは

$$
R_{\rm rms}=\frac{\sqrt{\left(99\right)}}{\left(100\right)}\simeq 0.099499
$$

である。

独立粒子の尺度からの補正は

$$
V_2/V_0=\frac{\left(0\right)}{\left(99\right)}=0
$$

である。

平均に対するゆらぎは一割程度になる。
:::

:::example[独立でも細かすぎる分解能]
小胞を小さくすると滑らかな密度が見えるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
N&=10000 \\p&=\frac{1}{10000} \\c&=0
\end{aligned}
$$

これらを順に代入する。

平均占有数は

$$
\mathbb E M=\left(10000\right) \left(\frac{1}{10000}\right)=1
$$

である。

対角分散は

$$
V_0=\left(10000\right) \left(\frac{1}{10000}\right) \left(1 - \left(\frac{1}{10000}\right)\right)=\frac{9999}{10000}
$$

である。

二体相関の分散寄与は

$$
V_2=\left(0\right) \left(10000\right) \left(\left(10000\right) - 1\right)=0
$$

である。

全占有数の分散は

$$
\operatorname{Var}M=\left(0\right) + \left(\frac{9999}{10000}\right)=\frac{9999}{10000}
$$

である。

密度の相対ゆらぎは

$$
R_{\rm rms}=\frac{\sqrt{\left(\frac{9999}{10000}\right)}}{\left(1\right)}\simeq 0.99995
$$

である。

独立粒子の尺度からの補正は

$$
V_2/V_0=\frac{\left(0\right)}{\left(\frac{9999}{10000}\right)}=0
$$

である。

平均個数は一個で、相対ゆらぎは一に近い。
:::

:::example[弱い正相関]
一対の相関が小さくても分散が増えるか求める。

入力値を次のように置く。

$$
\begin{aligned}
N&=10000 \\p&=\frac{1}{100} \\c&=\frac{1}{1000000}
\end{aligned}
$$

これらを順に代入する。

平均占有数は

$$
\mathbb E M=\left(10000\right) \left(\frac{1}{100}\right)=100
$$

である。

対角分散は

$$
V_0=\left(10000\right) \left(\frac{1}{100}\right) \left(1 - \left(\frac{1}{100}\right)\right)=99
$$

である。

二体相関の分散寄与は

$$
V_2=\left(10000\right) \left(\frac{1}{1000000}\right) \left(\left(10000\right) - 1\right)=\frac{9999}{100}
$$

である。

全占有数の分散は

$$
\operatorname{Var}M=\left(99\right) + \left(\frac{9999}{100}\right)=\frac{19899}{100}
$$

である。

密度の相対ゆらぎは

$$
R_{\rm rms}=\frac{\sqrt{\left(\frac{19899}{100}\right)}}{\left(100\right)}\simeq 0.14106
$$

である。

独立粒子の尺度からの補正は

$$
V_2/V_0=\frac{\left(\frac{9999}{100}\right)}{\left(99\right)}=\frac{101}{100}
$$

である。

対が多数あるため対角分散と同程度の寄与になる。
:::

:::example[ほぼ固定個数の小胞]
数の制約が生む負の相関を計算する。

入力値を次のように置く。

$$
\begin{aligned}
N&=10000 \\p&=\frac{1}{100} \\c&=- \frac{1}{1010000}
\end{aligned}
$$

これらを順に代入する。

平均占有数は

$$
\mathbb E M=\left(10000\right) \left(\frac{1}{100}\right)=100
$$

である。

対角分散は

$$
V_0=\left(10000\right) \left(\frac{1}{100}\right) \left(1 - \left(\frac{1}{100}\right)\right)=99
$$

である。

二体相関の分散寄与は

$$
V_2=\left(- \frac{1}{1010000}\right) \left(10000\right) \left(\left(10000\right) - 1\right)=-99
$$

である。

全占有数の分散は

$$
\operatorname{Var}M=\left(-99\right) + \left(99\right)=0
$$

である。

密度の相対ゆらぎは

$$
R_{\rm rms}=\frac{\sqrt{\left(0\right)}}{\left(100\right)}=0
$$

である。

独立粒子の尺度からの補正は

$$
V_2/V_0=\frac{\left(-99\right)}{\left(99\right)}=-1
$$

である。

総分散がゼロになる入力は、小胞に常に百個ある場合と一致する。
:::


## 線形化 Vlasov 方程式と Landau の初期値問題

一次元静電場、空間一様な平衡 $f_{s0}(v)=n_s F_s(v)$、$\int F_s\dd v=1$ とする。$k>0$ の Fourier 成分の摂動 $g_s(v,t)$ と $\phi(t)$ は
$$
\partial_tg_s+ikv g_s=ik\frac{q_sn_s}{m_s}\phi F_s',\qquad
\varepsilon_0k^2\phi=\sum_s q_s\int g_s\dd v
$$
を満たす。$\widetilde g(\omega)=\int_0^\infty e^{i\omega t}g(t)\dd t$、$\Im\omega$ を十分大きくとって境界項を消す。初期値を $g_{s0}$ とすると $\widetilde g_s=g_{s0}/[ik(v-\omega/k)]+(q_sn_s/m_s)\widetilde\phi F_s'/(v-\omega/k)$ である。従って
$$
\begin{aligned}
\varepsilon(\omega,k)&=1-\sum_s\frac{\omega_{ps}^2}{k^2}
\int_\R\frac{F_s'(v)}{v-\omega/k}\dd v,\\
\widetilde\phi(\omega)&=\frac{\sum_s q_s\int g_{s0}(v)/[i(kv-\omega)]\,\dd v}
{\varepsilon_0k^2\varepsilon(\omega,k)} .
\end{aligned}
\label{eq-plasma-kinetic-epsilon}
$$
逆変換は全ての成長極より上の水平線上で $\int e^{-i\omega t}\widetilde\phi\,\dd\omega/(2\pi)$ である。時間が正なら積分路を下へ下げる。応答の分母だけでなく、分子にも初期分布による速度の特異性がある。

:::example[位相混合]
力を除いた自由移流で $g_{s0}\propto e^{-v^2/(2v_t^2)}$ なら密度 Fourier 成分は初期値の $e^{-k^2v_t^2t^2/2}$ 倍になる。個々の $v$ の振幅は一定でも、速度間の位相差が密度の積分を小さくする。
:::

:::example[初期問題の実時間表示]
一次式を積分因子で解くと $g_s(t)=e^{-ikvt}g_{s0}+ik(q_sn_s/m_s)F_s'\int_0^t e^{-ikv(t-t')}\phi(t')\dd t'$ である。これを Poisson 方程式に入れると、時刻 $t$ の電位は過去の電位との畳み込みを含む。一つの分散根に初期値をすべて割り当てると、この自由輸送項が抜ける。
:::

### 初期分布から電場への Volterra 核

非磁化で一つの Maxwell 電子種を取り、平衡の平均流をゼロとする。空間波数 $k$ の摂動は $\partial_tf_1+ikvf_1=-(q/m)E_kF_0'$ に従う。積分因子を掛けて $f_1=e^{-ikvt}f_1(0)-(q/m)\int_0^t e^{-ikv(t-s)}E_k(s)F_0'\dd s$ を得る。Poisson 式 $ikE_k=q\int f_1\dd v/\varepsilon_0$ に入れ、$\int F_0'e^{-ikv\tau}\dd v=ik\tau n_0e^{-(kv_t\tau)^2/2}$ を使うと
$$
E_k(t)=E_{\rm free}(t)-\omega_p^2\int_0^t (t-s)
e^{-[kv_t(t-s)]^2/2} E_k(s)\dd s.
$$
冷たい極限では指数が一となるので二回微分すれば $\ddot E= \ddot E_{\rm free}-\omega_p^2E$ を得る。有限温度では履歴に重みが付き、二階の局所方程式に置き換えられない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $K$ | 波数 $kv_t/\omega_p$ | 1 |
| $s$ | 遅延時間 $\omega_p(t-s_{\rm past})$ | 1 |

熱速さで進む位相量を $X$ とする。波数と遅延時間の積である。

$$
X=K s
$$

記憶の Gaussian 因子を $F$ とする。自由飛行の速度積分から生じる。

$$
F=e^{- \frac{X^{2}}{2}}
$$

規格化した記憶核を $W$ とする。冷たい極限の線形な時間因子を掛ける。

$$
W=F s
$$

核が最大となる遅延を $s_{\rm peak}$ とする。$W'(s)=0$ を解く。

$$
s_{\rm peak}=\frac{1}{K}
$$

最大の記憶核を $W_{\rm peak}$ とする。$X=1$ で評価する。

$$
W_{\rm peak}=\frac{s_{\rm peak}}{e^{\frac{1}{2}}}
$$

核の規格化面積を $\int_0^\infty W\,\dd s$ とする。$X=k s$ と置換する。

$$
\int_0^\infty W\,\dd s=\frac{1}{K^{2}}
$$

:::example[長波長の短い遅延]
冷たい核とどれほど違うか求める。

入力値を次のように置く。

$$
\begin{aligned}
K&=\frac{1}{10} \\s&=1
\end{aligned}
$$

これらを順に代入する。

熱速さで進む位相量は

$$
X=\left(1\right) \left(\frac{1}{10}\right)=\frac{1}{10}
$$

である。

記憶の Gaussian 因子は

$$
F=e^{- \frac{\left(\frac{1}{10}\right)^{2}}{2}}\simeq 0.99501
$$

である。

規格化した記憶核は

$$
W=\left(1\right) \left(e^{- \frac{1}{200}}\right)\simeq 0.99501
$$

である。

核が最大となる遅延は

$$
s_{\rm peak}=\frac{1}{\left(\frac{1}{10}\right)}=10
$$

である。

最大の記憶核は

$$
W_{\rm peak}=\frac{\left(10\right)}{e^{\frac{1}{2}}}\simeq 6.0653
$$

である。

核の規格化面積は

$$
\int_0^\infty W\,\dd s=\frac{1}{\left(\frac{1}{10}\right)^{2}}=100
$$

である。

Gaussian 因子は一に近く、同じ時刻での冷たい極限に近い。
:::

:::example[長波長で長く待つ場合]
同じ波数でも冷模型が近くなくなるか算定する。

入力値を次のように置く。

$$
\begin{aligned}
K&=\frac{1}{10} \\s&=20
\end{aligned}
$$

これらを順に代入する。

熱速さで進む位相量は

$$
X=\left(20\right) \left(\frac{1}{10}\right)=2
$$

である。

記憶の Gaussian 因子は

$$
F=e^{- \frac{\left(2\right)^{2}}{2}}\simeq 0.13534
$$

である。

規格化した記憶核は

$$
W=\left(20\right) \left(e^{-2}\right)\simeq 2.7067
$$

である。

核が最大となる遅延は

$$
s_{\rm peak}=\frac{1}{\left(\frac{1}{10}\right)}=10
$$

である。

最大の記憶核は

$$
W_{\rm peak}=\frac{\left(10\right)}{e^{\frac{1}{2}}}\simeq 6.0653
$$

である。

核の規格化面積は

$$
\int_0^\infty W\,\dd s=\frac{1}{\left(\frac{1}{10}\right)^{2}}=100
$$

である。

遅延が十分長ければ位相混合で核が減る。
:::

:::example[熱速さで一波数長進む点]
履歴核のピークを確認する。

入力値を次のように置く。

$$
\begin{aligned}
K&=\frac{1}{2} \\s&=2
\end{aligned}
$$

これらを順に代入する。

熱速さで進む位相量は

$$
X=\left(2\right) \left(\frac{1}{2}\right)=1
$$

である。

記憶の Gaussian 因子は

$$
F=e^{- \frac{\left(1\right)^{2}}{2}}\simeq 0.60653
$$

である。

規格化した記憶核は

$$
W=\left(2\right) \left(e^{- \frac{1}{2}}\right)\simeq 1.2131
$$

である。

核が最大となる遅延は

$$
s_{\rm peak}=\frac{1}{\left(\frac{1}{2}\right)}=2
$$

である。

最大の記憶核は

$$
W_{\rm peak}=\frac{\left(2\right)}{e^{\frac{1}{2}}}\simeq 1.2131
$$

である。

核の規格化面積は

$$
\int_0^\infty W\,\dd s=\frac{1}{\left(\frac{1}{2}\right)^{2}}=4
$$

である。

この遅延で積の増加と指数の減少がつり合う。
:::

:::example[短波長の記憶]
遅い履歴と全積分の大きさを調べる。

入力値を次のように置く。

$$
\begin{aligned}
K&=2 \\s&=2
\end{aligned}
$$

これらを順に代入する。

熱速さで進む位相量は

$$
X=\left(2\right)^{2}=4
$$

である。

記憶の Gaussian 因子は

$$
F=e^{- \frac{\left(4\right)^{2}}{2}}\simeq 0.00033546
$$

である。

規格化した記憶核は

$$
W=\left(2\right) \left(e^{-8}\right)\simeq 0.00067093
$$

である。

核が最大となる遅延は

$$
s_{\rm peak}=\frac{1}{\left(2\right)}=\frac{1}{2}
$$

である。

最大の記憶核は

$$
W_{\rm peak}=\frac{\left(\frac{1}{2}\right)}{e^{\frac{1}{2}}}\simeq 0.30327
$$

である。

核の規格化面積は

$$
\int_0^\infty W\,\dd s=\frac{1}{\left(2\right)^{2}}=\frac{1}{4}
$$

である。

高波数では面積も減少するので、時間刻みは早いピークを分解する必要がある。
:::


## Landau 減衰

$\omega$ を実軸の上から近づけると分母は $v-\omega_r/k-i0$ となる。$1/(x-i\epsilon)=(x+i\epsilon)/(x^2+\epsilon^2)$ の虚部を積分すると $\pi\delta(x)$ が出る。実部は対称に特異点を除けば主値となる。従って解析接続では速度積分路を極の下に通し、
$$
\varepsilon_I(\omega_r,k)=-\pi\sum_s\frac{\omega_{ps}^2}{k^2}F_s'(\omega_r/k),\qquad
\gamma=-\frac{\varepsilon_I(\omega_r,k)}{\partial_{\omega_r}\varepsilon_R(\omega_r,k)}
\label{eq-plasma-landau-rate}
$$
を得る。第 2 式は単純な実根 $\varepsilon_R=0$ と $\abs\gamma\ll\abs{\omega_r}$ の近似である。平衡分布と初期摂動の解析性および無限遠での減衰により逆変換路を変形でき、減衰極の留数が非零なら、その極が指数減衰に寄与する。滑らかさだけでは極を越える変形は保証されない。:term[Landau 減衰]{en="Landau damping" yomi="らんだうげんすい"} は衝突がなくても起こる。電場のエネルギーは共鳴付近で粒子分布の変化へ移る。

:::example[符号の判定]
正の位相速度で $F_s'<0$、$\partial_\omega\varepsilon_R>0$ なら $\gamma<0$ である。分布に正の傾きを持つビームが加わると、その共鳴の寄与は成長側になる。
:::

### 実周波数の連続スペクトル

一種の規格化分布 $F$ を用い、初期値付きの式で時間依存を $e^{-iku t}$ としたときの一般化固有関数を考える。実数 $u$ に対し、分布の意味で
$$
g_u(v)=\frac{q n_0}{m}\phi_u
\operatorname{PV}\frac{F'(v)}{v-u}
+C_u\delta(v-u)
$$
と書ける。$(v-u)\operatorname{PV}[F'(v)/(v-u)]=F'(v)$ であり、$(v-u)\delta(v-u)=0$ なので線形 Vlasov 式を満たす。Poisson 式は $q C_u=\varepsilon_0 k^2\varepsilon_R(ku,k)\phi_u$ を与える。したがって実数 $u$ は、誘電関数の零点に限らず一般化固有関数を持つ。これを :term[van Kampen モード]{en="van Kampen mode" yomi="ばんかんぺんもーど"} と呼ぶ。

初期値が実周波数の重ね合わせで表され、電位の重み $A(u)$ が可積分なら、$\phi(t)=\int A(u)e^{-iku t}\dd u$ は $t\to\infty$ でゼロに近づく。速さが連続に分布するため、各成分の振幅が一定でも位相がそろわなくなる。Landau 処方は、同じ初期値を解析性が許す範囲で変形して減衰極の寄与を取り出す。連続スペクトルの完全性には関数空間と離散固有値の扱いが必要なので、ここでは実周波数の固有関数の構成までを示した。完全性を含む議論は :cite[nicholson-plasma] を参照する。

:::example[速度格子で位相が戻る時刻]
自由飛行を等間隔 $\Delta v$ の速度格子で計算すると、$t=2\pi/(|k|\Delta v)$ で隣の速度点の位相差が一周になる。有限和では全点の相対位相が初期値にそろう。この再帰時間は速度解像度に依存し、連続速度の物理的な減衰時間とは区別する。
:::

### 共鳴による電力の符号

実周波数 $\omega>0$ の静電波の電場と電流の振幅を $E,j$ とすると、$\sigma=-i\varepsilon_0\omega(\varepsilon-1)$ である。これは連続式の $j=(\omega/k)\rho_q$ と、誘起電荷からの誘電応答の定義を合わせて得る。従って実周波数で平均電力は $\frac12\varepsilon_0\omega \varepsilon_I\abs E^2$ であり、$\varepsilon_I>0$ の波では電場から粒子へエネルギーが渡る。位相速度の両側にまたがる共鳴粒子がいても、局所分布に傾きがあるので交換が相殺されない。

:::example[平坦な速度分布]
位相速度の近傍で $F_s$ が一定なら、その種の $F_s'(\omega_r/k)$ はゼロである。非共鳴の主値積分がゼロになる必要はないので、実周波数は依然として粒子によって変化する。
:::

### 共鳴粒子の密度と傾きの違い

平均 $U$、分散 $v_t^2$ の Maxwell 分布では $x=(v_\phi-U)/v_t$ として $F_0(v_\phi)/F_0(U)=e^{-x^2/2}$、$v_tF_0'(v_\phi)/F_0(U)=-xe^{-x^2/2}$ である。波のエネルギーが正で実根からのずれが小さいと、Landau 極の虚部の符号はこの速度勾配の符号に一致する。ただし誘電関数の周波数微分と波数を固定した比較である。$|x|=1$ で傾きの絶対値は最大になる。分布の中心では粒子数は最大であっても傾きが消える。この例で速度を固定して符号を比較することはできるが、流速変更により実根や波エネルギーも変わる系で成長率を一定の比例係数から予測しない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $x$ | 共鳴速度の中心からのずれを $v_t$ で割った値 | 1 |
| $K$ | 規格化波数 $k\lambda_D$ | 1 |

共鳴位置の分布の比を $F_r/F_c$ とする。正負のずれに対して対称である。

$$
F_r/F_c=e^{- \frac{x^{2}}{2}}
$$

規格化分布勾配を $G_r$ とする。電力交換には粒子数でなく勾配が入る。

$$
G_r=- \frac{F_r}{F_c} x
$$

規格化した分布の二階微分を $H_r$ とする。勾配の速度分解能を評価する。

$$
H_r=\frac{F_r}{F_c} \left(x^{2} - 1\right)
$$

最大傾きに対する符号付きの比を $G_r/\abs{G_{\max}}$ とする。$|x|=1$ を基準にする。

$$
G_r/\abs{G_{\max}}=G_{r} e^{\frac{1}{2}}
$$

誘電関数の虚部を $I_r$ とする。Maxwell 縦応答の $\sqrt\pi\,\zeta e^{-\zeta^2}/K^2$ である。

$$
I_r=- \frac{\sqrt{2} \sqrt{\pi} G_{r}}{2 K^{2}}
$$

変曲点での分布値に対する比を $F_r/F_{\rm inf}$ とする。分布の値と傾きを別々に比較する。

$$
F_r/F_{\rm inf}=\frac{F_r}{F_c} e^{\frac{1}{2}}
$$

:::example[分布の頂点]
共鳴に多数の粒子がいるときの虚部を求める。

入力値を次のように置く。

$$
\begin{aligned}
x&=0 \\K&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

共鳴位置の分布の比は

$$
F_r/F_c=e^{- \frac{\left(0\right)^{2}}{2}}=1
$$

である。

規格化分布勾配は

$$
G_r=- \left(0\right) \left(1\right)=0
$$

である。

規格化した分布の二階微分は

$$
H_r=\left(1\right) \left(\left(0\right)^{2} - 1\right)=-1
$$

である。

最大傾きに対する符号付きの比は

$$
G_r/\abs{G_{\max}}=\left(0\right) e^{\frac{1}{2}}=0
$$

である。

誘電関数の虚部は

$$
I_r=- \frac{\sqrt{2} \sqrt{\pi} \left(0\right)}{2 \left(\frac{1}{2}\right)^{2}}=0
$$

である。

変曲点での分布値に対する比は

$$
F_r/F_{\rm inf}=\left(1\right) e^{\frac{1}{2}}\simeq 1.6487
$$

である。

速度勾配がゼロなのでこの種からの一次吸収はゼロである。
:::

:::example[正の変曲点]
減衰を与える傾きの最大値を求める。

入力値を次のように置く。

$$
\begin{aligned}
x&=1 \\K&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

共鳴位置の分布の比は

$$
F_r/F_c=e^{- \frac{\left(1\right)^{2}}{2}}\simeq 0.60653
$$

である。

規格化分布勾配は

$$
G_r=- \left(1\right) \left(e^{- \frac{1}{2}}\right)\simeq -0.60653
$$

である。

規格化した分布の二階微分は

$$
H_r=\left(e^{- \frac{1}{2}}\right) \left(\left(1\right)^{2} - 1\right)=0
$$

である。

最大傾きに対する符号付きの比は

$$
G_r/\abs{G_{\max}}=\left(- \frac{1}{e^{\frac{1}{2}}}\right) e^{\frac{1}{2}}=-1
$$

である。

誘電関数の虚部は

$$
I_r=- \frac{\sqrt{2} \sqrt{\pi} \left(- \frac{1}{e^{\frac{1}{2}}}\right)}{2 \left(\frac{1}{2}\right)^{2}}\simeq 3.0407
$$

である。

変曲点での分布値に対する比は

$$
F_r/F_{\rm inf}=\left(e^{- \frac{1}{2}}\right) e^{\frac{1}{2}}=1
$$

である。

分布値は頂点より低いが傾きは最も急である。
:::

:::example[負の変曲点]
共鳴点が平均流より低い側の符号を求める。

入力値を次のように置く。

$$
\begin{aligned}
x&=-1 \\K&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

共鳴位置の分布の比は

$$
F_r/F_c=e^{- \frac{\left(-1\right)^{2}}{2}}\simeq 0.60653
$$

である。

規格化分布勾配は

$$
G_r=- \left(-1\right) \left(e^{- \frac{1}{2}}\right)\simeq 0.60653
$$

である。

規格化した分布の二階微分は

$$
H_r=\left(e^{- \frac{1}{2}}\right) \left(\left(-1\right)^{2} - 1\right)=0
$$

である。

最大傾きに対する符号付きの比は

$$
G_r/\abs{G_{\max}}=\left(e^{- \frac{1}{2}}\right) e^{\frac{1}{2}}=1
$$

である。

誘電関数の虚部は

$$
I_r=- \frac{\sqrt{2} \sqrt{\pi} \left(e^{- \frac{1}{2}}\right)}{2 \left(\frac{1}{2}\right)^{2}}\simeq -3.0407
$$

である。

変曲点での分布値に対する比は

$$
F_r/F_{\rm inf}=\left(e^{- \frac{1}{2}}\right) e^{\frac{1}{2}}=1
$$

である。

その種は正エネルギーの波へエネルギーを与え得る。
:::

:::example[高速の尾部]
遠い共鳴点で傾きが小さい理由を確認する。

入力値を次のように置く。

$$
\begin{aligned}
x&=4 \\K&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

共鳴位置の分布の比は

$$
F_r/F_c=e^{- \frac{\left(4\right)^{2}}{2}}\simeq 0.00033546
$$

である。

規格化分布勾配は

$$
G_r=- \left(4\right) \left(e^{-8}\right)\simeq -0.0013419
$$

である。

規格化した分布の二階微分は

$$
H_r=\left(e^{-8}\right) \left(\left(4\right)^{2} - 1\right)\simeq 0.0050319
$$

である。

最大傾きに対する符号付きの比は

$$
G_r/\abs{G_{\max}}=\left(- \frac{4}{e^{8}}\right) e^{\frac{1}{2}}\simeq -0.0022123
$$

である。

誘電関数の虚部は

$$
I_r=- \frac{\sqrt{2} \sqrt{\pi} \left(- \frac{4}{e^{8}}\right)}{2 \left(\frac{1}{2}\right)^{2}}\simeq 0.006727
$$

である。

変曲点での分布値に対する比は

$$
F_r/F_{\rm inf}=\left(e^{-8}\right) e^{\frac{1}{2}}\simeq 0.00055308
$$

である。

速度因子が増えても指数の減少が勝る。
:::


## プラズマ分散関数

$F(v)=(\sqrt{2\pi}v_t)^{-1}e^{-v^2/(2v_t^2)}$、$v_t=\sqrt{T/m}$ とする。:term[プラズマ分散関数]{en="plasma dispersion function" yomi="ぷらずまぶんさんかんすう"} は
$$
Z(\zeta)=\frac1{\sqrt\pi}\int_\R\frac{e^{-x^2}}{x-\zeta}\dd x,\qquad
\zeta=\frac{\omega}{\sqrt2 k v_t},\quad \Im\zeta>0
$$
を解析接続して定める。下半平面で積分路を極の下に保つには、実軸積分に $2i\sqrt\pi e^{-\zeta^2}$ を足す。これは実軸と、極を下から回り込む経路との差が反時計回りに極を囲むためである。上半平面の積分を微分後に部分積分すると $Z'=-2(1+\zeta Z)$、$Z(0)=i\sqrt\pi$ である。一次方程式を積分因子で解けば
$$
Z(\zeta)=e^{-\zeta^2}\left(i\sqrt\pi-2\int_0^\zeta e^{t^2}\dd t\right).
\label{eq-plasma-z-function}
$$
従って Maxwell 分布の応答は $\chi_s=(1+\zeta_s Z(\zeta_s))/(k^2\lambda_{Ds}^2)$ である。小 $\zeta$ では $Z=i\sqrt\pi-2\zeta+\order{\zeta^2}$。大きな実 $\zeta$ では分母を逆べきで展開してガウスモーメントを取る。実部は $-1/\zeta-1/(2\zeta^3)-3/(4\zeta^5)+\cdots$、虚部は $\sqrt\pi e^{-\zeta^2}$ である。

### Langmuir 波への適用

イオン固定、$K=k\lambda_{De}\ll1$ とする。実部の式から $\omega_r^2\simeq\omega_{pe}^2(1+3K^2)$、$\partial_\omega\varepsilon_R\simeq 2/\omega_{pe}$ である。$\varepsilon_I=\sqrt\pi\zeta e^{-\zeta^2}/K^2$ を使うと
$$
\frac{\gamma}{\omega_{pe}}\simeq-\sqrt{\frac{\pi}{8}}K^{-3}
\exp\left[-\frac{1}{2K^2}-\frac32\right].
$$
実周波数の $3K^2$ 補正が指数には有限の $-3/2$ を与える。現れる微小量をべき展開だけで求めることはできない。

:::example[原点の微分]
式 :eqref[eq-plasma-z-function] から $Z'(0)=-2$、もう一度微分して $Z''(0)=-2i\sqrt\pi$ である。
:::

:::example[実根の数値評価の手順]
$K=0.2$ の電子波では、Bohm–Gross 近似から $\omega_r^2/\omega_{pe}^2\simeq1.12$ なので $\zeta_r^2\simeq14$ である。弱減衰の先頭式では $-\sqrt{\pi/8}K^{-3}e^{-14}\simeq-6.5\times10^{-5}$ となる。減衰率の相対誤差が必要なら、先頭式の指数を使う代わりに式 :eqref[eq-plasma-z-function] の実部から $\varepsilon_R$ の根を数値的に求める。
:::

### プラズマ分散関数の実軸評価

$Z'(\zeta)=-2[1+\zeta Z(\zeta)]$ を使う。これは積分中で分母の微分を速度の全微分に移し、$\partial_x e^{-x^2}=-2xe^{-x^2}$ を代入して得る。$Z(0)=i\sqrt\pi$ なので、積分因子 $e^{\zeta^2}$ を掛けると実軸では $Z(x)=i\sqrt\pi e^{-x^2}-2D(x)$、$D(x)=e^{-x^2}\int_0^x e^{t^2}\dd t$ となる。誤差関数の定義で積分変数を純虚数に置くと $D(x)=\tfrac12\sqrt\pi e^{-x^2}\operatorname{erfi}(x)$ である。実部の大引数展開は $Z\simeq-1/x-1/(2x^3)-3/(4x^5)+\cdots$ であり、虚部の指数寄与は逆べき級数に現れない。以下の数値ではその虚部を独立に保持する。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $x$ | 非ゼロの実引数 | 1 |

Dawson 関数を $D(x)$ とする。原点からの実積分として求める。

$$
D(x)=\frac{\sqrt{\pi} e^{- x^{2}} \operatorname{erfi}{\left(x \right)}}{2}
$$

実部を $\Re Z$ とする。Dawson 関数の負の二倍である。

$$
\Re Z=- 2 D(x)
$$

虚部を $\Im Z$ とする。上側の解析接続で正になる。

$$
\Im Z=\sqrt{\pi} e^{- x^{2}}
$$

実部の逆べき近似を $Z_{\rm alg}$ とする。大引数の三項を残す。

$$
Z_{\rm alg}=- \frac{1}{x} - \frac{1}{2 x^{3}} - \frac{3}{4 x^{5}}
$$

三項近似の相対誤差を $\delta_R$ とする。実部のみを比較する。

$$
\delta_R=\left|{\frac{Z_{\rm alg} - \Re Z}{\Re Z}}\right|
$$

誘電応答に現れる実因子を $\Re(1+xZ)$ とする。大きな引数で二つの項が差し引かれる。

$$
\Re(1+xZ)=\Re Z x + 1
$$

:::example[原点に近い実引数]
大引数の近似を用いると何が起こるか確かめる。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

Dawson 関数は

$$
D(x)=\frac{\sqrt{\pi} e^{- \left(\frac{1}{2}\right)^{2}} \operatorname{erfi}{\left(\left(\frac{1}{2}\right) \right)}}{2}\simeq 0.42444
$$

である。

実部は

$$
\Re Z=- 2 \left(\frac{\sqrt{\pi} \operatorname{erfi}{\left(\frac{1}{2} \right)}}{2 e^{\frac{1}{4}}}\right)\simeq -0.84887
$$

である。

虚部は

$$
\Im Z=\sqrt{\pi} e^{- \left(\frac{1}{2}\right)^{2}}\simeq 1.3804
$$

である。

実部の逆べき近似は

$$
Z_{\rm alg}=- \frac{1}{\left(\frac{1}{2}\right)} - \frac{1}{2 \left(\frac{1}{2}\right)^{3}} - \frac{3}{4 \left(\frac{1}{2}\right)^{5}}=-30
$$

である。

三項近似の相対誤差は

$$
\delta_R=\left|{\frac{\left(- \frac{\sqrt{\pi} \operatorname{erfi}{\left(\frac{1}{2} \right)}}{e^{\frac{1}{4}}}\right) - \left(-30\right)}{\left(- \frac{\sqrt{\pi} \operatorname{erfi}{\left(\frac{1}{2} \right)}}{e^{\frac{1}{4}}}\right)}}\right|\simeq 34.341
$$

である。

誘電応答に現れる実因子は

$$
\Re(1+xZ)=\left(- \frac{\sqrt{\pi} \operatorname{erfi}{\left(\frac{1}{2} \right)}}{e^{\frac{1}{4}}}\right) \left(\frac{1}{2}\right) + 1\simeq 0.57556
$$

である。

逆べき近似の誤差が大きいので実積分を使う。
:::

:::example[中間の実引数]
実部と虚部を同時に求める。

入力値を次のように置く。

$$
\begin{aligned}
x&=2
\end{aligned}
$$

これらを順に代入する。

Dawson 関数は

$$
D(x)=\frac{\sqrt{\pi} e^{- \left(2\right)^{2}} \operatorname{erfi}{\left(\left(2\right) \right)}}{2}\simeq 0.30134
$$

である。

実部は

$$
\Re Z=- 2 \left(\frac{\sqrt{\pi} \operatorname{erfi}{\left(2 \right)}}{2 e^{4}}\right)\simeq -0.60268
$$

である。

虚部は

$$
\Im Z=\sqrt{\pi} e^{- \left(2\right)^{2}}\simeq 0.032464
$$

である。

実部の逆べき近似は

$$
Z_{\rm alg}=- \frac{1}{\left(2\right)} - \frac{1}{2 \left(2\right)^{3}} - \frac{3}{4 \left(2\right)^{5}}=- \frac{75}{128}
$$

である。

三項近似の相対誤差は

$$
\delta_R=\left|{\frac{\left(- \frac{75}{128}\right) - \left(- \frac{\sqrt{\pi} \operatorname{erfi}{\left(2 \right)}}{e^{4}}\right)}{\left(- \frac{\sqrt{\pi} \operatorname{erfi}{\left(2 \right)}}{e^{4}}\right)}}\right|\simeq 0.027781
$$

である。

誘電応答に現れる実因子は

$$
\Re(1+xZ)=\left(- \frac{\sqrt{\pi} \operatorname{erfi}{\left(2 \right)}}{e^{4}}\right) \left(2\right) + 1\simeq -0.20536
$$

である。

虚部は小さくなってもゼロではない。
:::

:::example[高速相速度の極限]
実部の三項展開と正確な値を比較する。

入力値を次のように置く。

$$
\begin{aligned}
x&=5
\end{aligned}
$$

これらを順に代入する。

Dawson 関数は

$$
D(x)=\frac{\sqrt{\pi} e^{- \left(5\right)^{2}} \operatorname{erfi}{\left(\left(5\right) \right)}}{2}\simeq 0.10213
$$

である。

実部は

$$
\Re Z=- 2 \left(\frac{\sqrt{\pi} \operatorname{erfi}{\left(5 \right)}}{2 e^{25}}\right)\simeq -0.20427
$$

である。

虚部は

$$
\Im Z=\sqrt{\pi} e^{- \left(5\right)^{2}}\simeq 2.4616\times 10^{-11}
$$

である。

実部の逆べき近似は

$$
Z_{\rm alg}=- \frac{1}{\left(5\right)} - \frac{1}{2 \left(5\right)^{3}} - \frac{3}{4 \left(5\right)^{5}}=- \frac{2553}{12500}
$$

である。

三項近似の相対誤差は

$$
\delta_R=\left|{\frac{\left(- \frac{2553}{12500}\right) - \left(- \frac{\sqrt{\pi} \operatorname{erfi}{\left(5 \right)}}{e^{25}}\right)}{\left(- \frac{\sqrt{\pi} \operatorname{erfi}{\left(5 \right)}}{e^{25}}\right)}}\right|\simeq 0.0001378
$$

である。

誘電応答に現れる実因子は

$$
\Re(1+xZ)=\left(- \frac{\sqrt{\pi} \operatorname{erfi}{\left(5 \right)}}{e^{25}}\right) \left(5\right) + 1\simeq -0.021341
$$

である。

実部を近似できても虚部は別の指数で評価する。
:::

:::example[負の引数]
実部と虚部の対称性を確かめる。

入力値を次のように置く。

$$
\begin{aligned}
x&=-2
\end{aligned}
$$

これらを順に代入する。

Dawson 関数は

$$
D(x)=\frac{\sqrt{\pi} e^{- \left(-2\right)^{2}} \operatorname{erfi}{\left(\left(-2\right) \right)}}{2}\simeq -0.30134
$$

である。

実部は

$$
\Re Z=- 2 \left(- \frac{\sqrt{\pi} \operatorname{erfi}{\left(2 \right)}}{2 e^{4}}\right)\simeq 0.60268
$$

である。

虚部は

$$
\Im Z=\sqrt{\pi} e^{- \left(-2\right)^{2}}\simeq 0.032464
$$

である。

実部の逆べき近似は

$$
Z_{\rm alg}=- \frac{1}{\left(-2\right)} - \frac{1}{2 \left(-2\right)^{3}} - \frac{3}{4 \left(-2\right)^{5}}=\frac{75}{128}
$$

である。

三項近似の相対誤差は

$$
\delta_R=\left|{\frac{\left(\frac{75}{128}\right) - \left(\frac{\sqrt{\pi} \operatorname{erfi}{\left(2 \right)}}{e^{4}}\right)}{\left(\frac{\sqrt{\pi} \operatorname{erfi}{\left(2 \right)}}{e^{4}}\right)}}\right|\simeq 0.027781
$$

である。

誘電応答に現れる実因子は

$$
\Re(1+xZ)=\left(-2\right) \left(\frac{\sqrt{\pi} \operatorname{erfi}{\left(2 \right)}}{e^{4}}\right) + 1\simeq -0.20536
$$

である。

実部は奇関数、虚部は偶関数である。
:::


## イオン音波の Landau 減衰

$T_e\gg T_i$、$k\lambda_{De}\ll1$ で正のイオン音波を考える。$\zeta_e\ll1$、$\zeta_i\gg1$ なので、実部の主要項は $1+1/(k^2\lambda_{De}^2)-\omega_{pi}^2/\omega^2$。微分は $2\omega_{pi}^2/\omega_r^3$ である。虚部に両種を入れて
$$
\frac{\gamma}{\omega_r}\simeq
-\sqrt{\frac{\pi}{8}}
\left[\sqrt{\frac{m_e}{m_i}}+
\left(\frac{T_e}{T_i}\right)^{3/2}e^{-\omega_r^2/(2k^2v_{ti}^2)}\right]
$$
を得る。ここでは虚部の前置係数を $\omega_r/k\simeq\sqrt{T_e/m_i}$ で簡約した。指数には高次の周波数補正が影響するので、定量計算では実誘電関数の根を使う。$T_e/T_i$ が小さい場合には弱減衰近似そのものを見直す。

:::example[電子だけが効く極限]
$T_i/T_e\to0$ でイオン項は指数的に消える。陽子なら $\abs{\gamma}/\omega_r\simeq\sqrt{\pi/8}/\sqrt{1836}\simeq0.0146$ となる。
:::

:::example[有限イオン温度の項を比べる]
$\omega_r/k$ を冷イオン音速で置き、$T_e/T_i=30$ とする。イオン項の前置係数と指数を掛けると $30^{3/2}e^{-15}\simeq5.0\times10^{-5}$ である。陽子での電子項 $1/\sqrt{1836}\simeq2.3\times10^{-2}$ より小さい。この比較は冷イオン音速を使った目安で、イオン温度による実周波数変化を入れるとイオン項の指数も変わる。
:::

### イオン音波で二種の共鳴を比べる

一価水素の等方な Maxwell 分布で平均流はゼロとする。$\tau=T_i/T_e$、$K=k\lambda_{De}$ と置く。電子は Boltzmann、イオンの圧力補正を一次元断熱とすれば $v_\phi^2=(T_e/m_i)[(1+K^2)^{-1}+3\tau]$ である。この位相速度を運動論の弱減衰推定に使う。$\zeta_i^2=v_\phi^2/(2T_i/m_i)$、$\zeta_e^2=v_\phi^2/(2T_e/m_e)$ である。誘電関数の種ごとの虚部は $\sqrt\pi\,\zeta_s e^{-\zeta_s^2}/(k\lambda_{Ds})^2$ だから、その比には $\lambda_{De}^2/\lambda_{Di}^2=1/\tau$ が入る。二つの虚部の比較は実周波数が与えられているときの診断であり、大きな減衰の場合には複素零点を取り直す。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $\tau$ | 温度比 $T_i/T_e$ | 1 |
| $K$ | $k\lambda_{De}$ | 1 |
| $r$ | 質量比 $m_e/m_i$ | 1 |

電子温度単位での位相速度二乗を $W$ とする。遮蔽とイオン圧の双方を含む近似である。

$$
W=3 \tau + \frac{1}{K^{2} + 1}
$$

イオンの共鳴引数を $\zeta_i$ とする。イオン熱速さで規格化する。

$$
\zeta_i=\frac{\sqrt{2} \sqrt{\frac{W}{\tau}}}{2}
$$

電子の共鳴引数を $\zeta_e$ とする。電子の軽さで引数が小さくなる。

$$
\zeta_e=\frac{\sqrt{2} \sqrt{W r}}{2}
$$

イオンの虚部寄与を $I_i$ とする。イオン Debye 長を使う。

$$
I_i=\frac{\sqrt{\pi} \zeta_{i} e^{- \zeta_{i}^{2}}}{K^{2} \tau}
$$

電子の虚部寄与を $I_e$ とする。電子は分布の中心付近で共鳴する。

$$
I_e=\frac{\sqrt{\pi} \zeta_{e} e^{- \zeta_{e}^{2}}}{K^{2}}
$$

二つの虚部の比を $I_i/I_e$ とする。同じ実周波数で比較する。

$$
I_i/I_e=\frac{I_{i}}{I_{e}}
$$

:::example[冷たいイオン]
イオンの共鳴が指数的に小さいか確かめる。

入力値を次のように置く。

$$
\begin{aligned}
\tau&=\frac{1}{100} \\K&=\frac{1}{10} \\r&=\frac{1}{1836}
\end{aligned}
$$

これらを順に代入する。

電子温度単位での位相速度二乗は

$$
W=3 \left(\frac{1}{100}\right) + \frac{1}{\left(\frac{1}{10}\right)^{2} + 1}=\frac{10303}{10100}
$$

である。

イオンの共鳴引数は

$$
\zeta_i=\frac{\sqrt{2} \sqrt{\frac{\left(\frac{10303}{10100}\right)}{\left(\frac{1}{100}\right)}}}{2}\simeq 7.1418
$$

である。

電子の共鳴引数は

$$
\zeta_e=\frac{\sqrt{2} \sqrt{\left(\frac{10303}{10100}\right) \left(\frac{1}{1836}\right)}}{2}\simeq 0.016667
$$

である。

イオンの虚部寄与は

$$
I_i=\frac{\sqrt{\pi} \left(\frac{\sqrt{2081206}}{202}\right) e^{- \left(\frac{\sqrt{2081206}}{202}\right)^{2}}}{\left(\frac{1}{100}\right) \left(\frac{1}{10}\right)^{2}}\simeq 8.9374\times 10^{-18}
$$

である。

電子の虚部寄与は

$$
I_e=\frac{\sqrt{\pi} \left(\frac{\sqrt{106141506}}{618120}\right) e^{- \left(\frac{\sqrt{106141506}}{618120}\right)^{2}}}{\left(\frac{1}{10}\right)^{2}}\simeq 2.9534
$$

である。

二つの虚部の比は

$$
I_i/I_e=\frac{\left(\frac{5000 \sqrt{2081206} \sqrt{\pi}}{101 e^{\frac{10303}{202}}}\right)}{\left(\frac{5 \sqrt{106141506} \sqrt{\pi}}{30906 e^{\frac{10303}{37087200}}}\right)}\simeq 3.0261\times 10^{-18}
$$

である。

イオン圧の補正を加えても位相速度はイオン熱速さより大きい。
:::

:::example[温度差が小さくなる場合]
電子とイオンでどちらの寄与が大きいかを求める。

入力値を次のように置く。

$$
\begin{aligned}
\tau&=\frac{1}{5} \\K&=\frac{1}{10} \\r&=\frac{1}{1836}
\end{aligned}
$$

これらを順に代入する。

電子温度単位での位相速度二乗は

$$
W=3 \left(\frac{1}{5}\right) + \frac{1}{\left(\frac{1}{10}\right)^{2} + 1}=\frac{803}{505}
$$

である。

イオンの共鳴引数は

$$
\zeta_i=\frac{\sqrt{2} \sqrt{\frac{\left(\frac{803}{505}\right)}{\left(\frac{1}{5}\right)}}}{2}\simeq 1.9938
$$

である。

電子の共鳴引数は

$$
\zeta_e=\frac{\sqrt{2} \sqrt{\left(\frac{1}{1836}\right) \left(\frac{803}{505}\right)}}{2}\simeq 0.020809
$$

である。

イオンの虚部寄与は

$$
I_i=\frac{\sqrt{\pi} \left(\frac{\sqrt{162206}}{202}\right) e^{- \left(\frac{\sqrt{162206}}{202}\right)^{2}}}{\left(\frac{1}{10}\right)^{2} \left(\frac{1}{5}\right)}\simeq 33.174
$$

である。

電子の虚部寄与は

$$
I_e=\frac{\sqrt{\pi} \left(\frac{\sqrt{41362530}}{309060}\right) e^{- \left(\frac{\sqrt{41362530}}{309060}\right)^{2}}}{\left(\frac{1}{10}\right)^{2}}\simeq 3.6868
$$

である。

二つの虚部の比は

$$
I_i/I_e=\frac{\left(\frac{250 \sqrt{162206} \sqrt{\pi}}{101 e^{\frac{803}{202}}}\right)}{\left(\frac{5 \sqrt{41362530} \sqrt{\pi}}{15453 e^{\frac{803}{1854360}}}\right)}\simeq 8.9981
$$

である。

イオン共鳴が競合し、冷イオンに基づく弱減衰の近似を再検討する。
:::

:::example[Debye 遮蔽が強い波数]
同じ温度でも波数によって共鳴速度が変わることを調べる。

入力値を次のように置く。

$$
\begin{aligned}
\tau&=\frac{1}{10} \\K&=1 \\r&=\frac{1}{1836}
\end{aligned}
$$

これらを順に代入する。

電子温度単位での位相速度二乗は

$$
W=3 \left(\frac{1}{10}\right) + \frac{1}{\left(1\right)^{2} + 1}=\frac{4}{5}
$$

である。

イオンの共鳴引数は

$$
\zeta_i=\frac{\sqrt{2} \sqrt{\frac{\left(\frac{4}{5}\right)}{\left(\frac{1}{10}\right)}}}{2}=2
$$

である。

電子の共鳴引数は

$$
\zeta_e=\frac{\sqrt{2} \sqrt{\left(\frac{1}{1836}\right) \left(\frac{4}{5}\right)}}{2}\simeq 0.01476
$$

である。

イオンの虚部寄与は

$$
I_i=\frac{\sqrt{\pi} \left(2\right) e^{- \left(2\right)^{2}}}{\left(1\right)^{2} \left(\frac{1}{10}\right)}\simeq 0.64927
$$

である。

電子の虚部寄与は

$$
I_e=\frac{\sqrt{\pi} \left(\frac{\sqrt{510}}{1530}\right) e^{- \left(\frac{\sqrt{510}}{1530}\right)^{2}}}{\left(1\right)^{2}}\simeq 0.026156
$$

である。

二つの虚部の比は

$$
I_i/I_e=\frac{\left(\frac{20 \sqrt{\pi}}{e^{4}}\right)}{\left(\frac{\sqrt{510} \sqrt{\pi}}{1530 e^{\frac{1}{4590}}}\right)}\simeq 24.823
$$

である。

位相速度が下がり、イオンの指数的抑制が弱まる。
:::

:::example[等質量で異なる温度]
電子 Boltzmann の仮定が成り立つか形式値で点検する。

入力値を次のように置く。

$$
\begin{aligned}
\tau&=\frac{1}{10} \\K&=\frac{1}{10} \\r&=1
\end{aligned}
$$

これらを順に代入する。

電子温度単位での位相速度二乗は

$$
W=3 \left(\frac{1}{10}\right) + \frac{1}{\left(\frac{1}{10}\right)^{2} + 1}=\frac{1303}{1010}
$$

である。

イオンの共鳴引数は

$$
\zeta_i=\frac{\sqrt{2} \sqrt{\frac{\left(\frac{1303}{1010}\right)}{\left(\frac{1}{10}\right)}}}{2}\simeq 2.5398
$$

である。

電子の共鳴引数は

$$
\zeta_e=\frac{\sqrt{2} \sqrt{\left(1\right) \left(\frac{1303}{1010}\right)}}{2}\simeq 0.80315
$$

である。

イオンの虚部寄与は

$$
I_i=\frac{\sqrt{\pi} \left(\frac{\sqrt{263206}}{202}\right) e^{- \left(\frac{\sqrt{263206}}{202}\right)^{2}}}{\left(\frac{1}{10}\right)^{3}}\simeq 7.1114
$$

である。

電子の虚部寄与は

$$
I_e=\frac{\sqrt{\pi} \left(\frac{\sqrt{658015}}{1010}\right) e^{- \left(\frac{\sqrt{658015}}{1010}\right)^{2}}}{\left(\frac{1}{10}\right)^{2}}\simeq 74.684
$$

である。

二つの虚部の比は

$$
I_i/I_e=\frac{\left(\frac{500 \sqrt{263206} \sqrt{\pi}}{101 e^{\frac{1303}{202}}}\right)}{\left(\frac{10 \sqrt{658015} \sqrt{\pi}}{101 e^{\frac{1303}{2020}}}\right)}\simeq 0.09522
$$

である。

電子引数は一より十分小さくなく、この位相速度を使ったイオン音波の推定は適用できない。
:::


## Bernstein–Greene–Kruskal モード

電位が波とともに進むとき、その座標 $\xi=x-Ut$ に移れば静的な問題になる。粒子のエネルギー $\mathcal H_s=m_s(v-U)^2/2+q_s\phi(\xi)$ は軌道に沿って一定である。従って $f_s=G_s(\mathcal H_s,\sigma)$ は定常の Vlasov 方程式を満たす。離れた通過軌道には符号 $\sigma=\sgn(v-U)$ を付け、捕捉軌道には独立な $G_s$ を与えられる。Poisson 方程式も満たすものが :term[Bernstein–Greene–Kruskal モード]{en="Bernstein–Greene–Kruskal mode" yomi="ばーんすたいんぐりーんくるすかるもーど"} である。

密度が $\phi$ の一価関数になる分布を選ぶと、$\phi''=-\rho_q(\phi)/\varepsilon_0$ である。$\phi'$ を掛けて積分すれば $(\phi')^2/2+\int^\phi \rho_q(y)\dd y/\varepsilon_0=C$ となる。周期波はこの有効ポテンシャル内の往復に対応する。逆に与えた電位から正の捕捉分布を再構成するには、通過粒子の密度を引いた Abel 積分を反転し、分離軌道での整合も確かめる。任意の電位が物理的な分布を持つことは保証されない。その構成の詳論は :cite[nicholson-plasma] を参照する。

:::example[捕捉域の幅]
電子の最大電位と分離軌道の電位の差を $\Delta\phi>0$ とすると、最大電位の点で $\abs{v-U}<\sqrt{2e\Delta\phi/m_e}$ が捕捉域である。電位差を四倍にすると速度幅は二倍になる。
:::

:::example[通過軌道と捕捉軌道]
ポテンシャルが周期的に変わる場合、山より高いエネルギーの通過粒子は相対速度の符号を変えない。それぞれの符号に異なる定常分布を置ける。井戸に捕まった軌道は速度ゼロの反射点を通って符号を変えるため、同じ閉軌道の二つの枝には同じ分布を与える。
:::

### Abel 反転による構成

具体的に、進行座標で規格化したポテンシャル $\phi\ge0$ と電子の相対エネルギー $E=v^2/2-\phi$ を使う。長さ、速度、密度を選び Poisson 方程式を $\phi''=n_e-n_i$ の形にし、一定のイオン背景を考える。通過粒子は $E>0$ の左右で同じ分布を置き、捕捉粒子は $E=-u$、$0<u<\phi$ の枝に $G_t(u)$ を置く。$v=\pm\sqrt{2(\phi-u)}$ の両側を数えると捕捉密度は
$$
n_t(\phi)=\sqrt2\int_0^\phi\frac{G_t(u)}{\sqrt{\phi-u}}\dd u.
$$
これに $(s-\phi)^{-1/2}$ を掛けて $\phi=0$ から $s$ まで積分する。順序を交換すると、内部積分は
$$
\int_u^s\frac{\dd\phi}{\sqrt{(\phi-u)(s-\phi)}}=\pi
$$
となる。$\phi=u+(s-u)\sin^2\theta$ と置けば $2\int_0^{\pi/2}\dd\theta=\pi$ である。従って
$$
G_t(s)=\frac{1}{\sqrt2\pi}\dv{}{s}\int_0^s
\frac{n_t(\phi)}{\sqrt{s-\phi}}\dd\phi .
$$
$n_t(0)=0$ で $n_t$ が微分可能なら、積分を $s-\phi$ で書き直して $G_t(s)=(\sqrt2\pi)^{-1}\int_0^s n_t'(\phi)/\sqrt{s-\phi}\dd\phi$ とできる。従って $n_t'\ge0$ は正の捕捉分布を得るための十分条件である。

:::example[周期電位の自己無撞着な例]
$\phi(\xi)=\phi_m[1+\cos(k\xi)]/2$ とする。要求される電子密度は $n_e(\phi)=n_i+k^2(\phi_m/2-\phi)$ である。通過粒子は規格化した単一エネルギー $E_0>0$ に集中させ、$\phi=0$ での密度を $N=n_i+k^2\phi_m/2$ とする。各方向の定常分布を $N\sqrt{E_0/2}\delta(E-E_0)$ と置くと、速度積分による密度は $n_p(\phi)=N/\sqrt{1+\phi/E_0}$ になる。必要な捕捉密度
$$
n_t(\phi)=N\left(1-\frac1{\sqrt{1+\phi/E_0}}\right)-k^2\phi
$$
は原点でゼロである。$k^2\le N[2E_0(1+\phi_m/E_0)^{3/2}]^{-1}$ なら $0\le\phi\le\phi_m$ で $n_t'\ge0$ であり、反転で非負の捕捉分布が得られる。例えば $n_i=1$、$\phi_m=E_0=1$、$k^2=1/10$ なら右辺は $21/(80\sqrt2)>1/10$ である。電荷密度の和を Poisson 方程式に戻せば指定した余弦電位を再現する。通過分布のデルタは、軌道上に支えを持つ弱い Vlasov 解として理解する。
:::

### 捕捉粒子の密度を Abel 積分で組む

進行座標で $H=mv^2/2-\psi$ とし、$\psi\ge0$ をポテンシャル井戸の深さとする。捕捉分布 $F_t(H)$ は $-\psi\le H<0$ にある。正負の速度を足し、$\dd v=\dd H/\sqrt{2m(H+\psi)}$ と置換して $n_t(\psi)=\sqrt{2/m}\int_{-\psi}^0 F_t(H)(H+\psi)^{-1/2}\dd H$ を得る。$F_t(H)=A+B(-H)$ なら $n_t=\sqrt{2/m}[2A\sqrt\psi+(4/3)B\psi^{3/2}]$ である。逆に一つの密度の振幅が分かっても $F_t$ の形全体は分からない。Poisson 方程式の右辺の残りを通過密度で補い、全粒子密度が非負であるかを確認する。ここでは指定の捕捉分布から密度を得る計算を行う。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $x$ | 井戸深さを一定エネルギーで割った値 | 1 |
| $A$ | 捕捉分布の定数項の規格化値 | 1 |
| $B$ | 捕捉分布の一次項の規格化値 | 1 |

捕捉速度の上限を $V_t$ とする。速度単位を $\sqrt{2\psi_*/m}$ とする。

$$
V_t=\sqrt{x}
$$

定数分布からの密度寄与を $N_A$ とする。密度の単位は上の速度単位と分布単位の積である。

$$
N_A=2 A V_{t}
$$

一次分布からの密度寄与を $N_B$ とする。Abel 積分の異なるべきで増える。

$$
N_B=\frac{4 B V_{t} x}{3}
$$

捕捉密度を $N_t$ とする。二つの積分を加える。

$$
N_t=N_{A} + N_{B}
$$

井戸の底での捕捉分布を $F_t(-x)$ とする。非負性はこの値と $A$ で確認する。

$$
F_t(-x)=A + B x
$$

一次項の密度寄与の比を $N_B/N_t$ とする。負なら平坦分布に対する欠損を表す。

$$
N_B/N_t=\frac{N_{B}}{N_{t}}
$$

:::example[平坦な捕捉分布]
井戸の深さが小さい領域の密度を求める。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{4} \\A&=1 \\B&=0
\end{aligned}
$$

これらを順に代入する。

捕捉速度の上限は

$$
V_t=\sqrt{\left(\frac{1}{4}\right)}=\frac{1}{2}
$$

である。

定数分布からの密度寄与は

$$
N_A=2 \left(1\right) \left(\frac{1}{2}\right)=1
$$

である。

一次分布からの密度寄与は

$$
N_B=\frac{4 \left(0\right) \left(\frac{1}{2}\right) \left(\frac{1}{4}\right)}{3}=0
$$

である。

捕捉密度は

$$
N_t=\left(0\right) + \left(1\right)=1
$$

である。

井戸の底での捕捉分布は

$$
F_t(-x)=\left(0\right) \left(\frac{1}{4}\right) + \left(1\right)=1
$$

である。

一次項の密度寄与の比は

$$
N_B/N_t=\frac{\left(0\right)}{\left(1\right)}=0
$$

である。

捕捉密度は速度幅に比例する。
:::

:::example[内側に粒子を多く置く場合]
同じ分離軌道で内部を増やしたときの密度を求める。

入力値を次のように置く。

$$
\begin{aligned}
x&=1 \\A&=1 \\B&=1
\end{aligned}
$$

これらを順に代入する。

捕捉速度の上限は

$$
V_t=\sqrt{\left(1\right)}=1
$$

である。

定数分布からの密度寄与は

$$
N_A=2 \left(1\right)^{2}=2
$$

である。

一次分布からの密度寄与は

$$
N_B=\frac{4 \left(1\right)^{3}}{3}=\frac{4}{3}
$$

である。

捕捉密度は

$$
N_t=\left(2\right) + \left(\frac{4}{3}\right)=\frac{10}{3}
$$

である。

井戸の底での捕捉分布は

$$
F_t(-x)=\left(1\right)^{2} + \left(1\right)=2
$$

である。

一次項の密度寄与の比は

$$
N_B/N_t=\frac{\left(\frac{4}{3}\right)}{\left(\frac{10}{3}\right)}=\frac{2}{5}
$$

である。

一次項は深い井戸ほど相対的に大きい。
:::

:::example[捕捉粒子の欠損]
分布の非負性を保って密度が減る例を計算する。

入力値を次のように置く。

$$
\begin{aligned}
x&=1 \\A&=1 \\B&=- \frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

捕捉速度の上限は

$$
V_t=\sqrt{\left(1\right)}=1
$$

である。

定数分布からの密度寄与は

$$
N_A=2 \left(1\right)^{2}=2
$$

である。

一次分布からの密度寄与は

$$
N_B=\frac{4 \left(- \frac{1}{2}\right) \left(1\right)^{2}}{3}=- \frac{2}{3}
$$

である。

捕捉密度は

$$
N_t=\left(- \frac{2}{3}\right) + \left(2\right)=\frac{4}{3}
$$

である。

井戸の底での捕捉分布は

$$
F_t(-x)=\left(- \frac{1}{2}\right) \left(1\right) + \left(1\right)=\frac{1}{2}
$$

である。

一次項の密度寄与の比は

$$
N_B/N_t=\frac{\left(- \frac{2}{3}\right)}{\left(\frac{4}{3}\right)}=- \frac{1}{2}
$$

である。

欠損は負の分布そのものを意味しない。
:::

:::example[しきい値で消える捕捉種]
分離軌道でゼロの場合の深さ依存を調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{4} \\A&=0 \\B&=1
\end{aligned}
$$

これらを順に代入する。

捕捉速度の上限は

$$
V_t=\sqrt{\left(\frac{1}{4}\right)}=\frac{1}{2}
$$

である。

定数分布からの密度寄与は

$$
N_A=2 \left(0\right) \left(\frac{1}{2}\right)=0
$$

である。

一次分布からの密度寄与は

$$
N_B=\frac{4 \left(1\right) \left(\frac{1}{2}\right) \left(\frac{1}{4}\right)}{3}=\frac{1}{6}
$$

である。

捕捉密度は

$$
N_t=\left(0\right) + \left(\frac{1}{6}\right)=\frac{1}{6}
$$

である。

井戸の底での捕捉分布は

$$
F_t(-x)=\left(0\right) + \left(1\right) \left(\frac{1}{4}\right)=\frac{1}{4}
$$

である。

一次項の密度寄与の比は

$$
N_B/N_t=1
$$

である。

定数項がないので先頭は $\psi^{3/2}$ になる。
:::


## 演習問題

::::exercise{level=basic}
場のない系で初期値が $g(x,v)$ の解を書け。

:::solution
特性は $x(t)-vt=x(0)$、$v$ 一定なので解は $g(x-vt,v)$ である。
:::
::::

::::exercise{level=basic}
$F$ が偶関数なら $F'$ はどの対称性をもつか。

:::solution
$F(-v)=F(v)$ を微分すると $F'(-v)=-F'(v)$ である。従って微分は奇関数である。
:::
::::

::::exercise{level=basic}
$Z'(0)$ を求めよ。

:::solution
$Z'=-2(1+\zeta Z)$ に $\zeta=0$ を入れて $-2$ となる。
:::
::::

::::exercise{level=basic}
$k\lambda_{De}\ll1$、$T_i/T_e\to0$ の音波で $m_i/m_e=2500$ のとき $\gamma/\omega_r$ の見積もりを求めよ。

:::solution
電子項だけ残るので $-\sqrt{\pi/8}/50$ である。
:::
::::

::::exercise{level=advanced}
平均ゼロ、分散 $v_t^2$ の Maxwell 分布 $F$ に対して $\int_{\R}v^3 F'(v)\dd v$ を求めよ。

:::solution
部分積分により $[v^3F]_{-\infty}^{\infty}-3\int v^2F\dd v=-3v_t^2$ である。
:::
::::

::::exercise{level=advanced}
$\varepsilon_I=0.02$、$\partial_\omega\varepsilon_R=0.01\,\mathrm s$ の弱減衰根の $\gamma$ を求めよ。

:::solution
一次補正式 :eqref[eq-plasma-landau-rate] により $\gamma=-2\,\mathrm{s^{-1}}$ である。これは $\abs\gamma$ が $\omega_r$ より小さい場合の値である。
:::
::::

::::exercise{level=advanced}
分散関数の第三微分を原点で求めよ。

:::solution
$Z''=-2Z-2\zeta Z'$、$Z'''=-4Z'-2\zeta Z''$ なので $Z'''(0)=8$ である。
:::
::::

::::exercise{level=advanced}
$\phi(\xi)$ が進行波であるとき $\mathcal H_s$ の保存を証明せよ。

:::solution
$\dot\xi=v-U$、$m_s\dot v=-q_s\phi'$ なので $\dot{\mathcal H}_s=m_s(v-U)\dot v+q_s\phi'(v-U)=0$ である。
:::
::::

## 確認テスト

::::quiz
単一の減衰電場の位相混合に伴う Vlasov 分布について正しいものはどれか。

- [x] 分布は位相空間で細かな構造を作りうる。
- [ ] 衝突項が必ず非ゼロになる。
- [ ] 各軌道上で分布が指数減衰する。

:::explanation
式 :eqref[eq-plasma-vlasov] により軌道上で $f_s$ は一定である。異なる速度の位相がずれ、積分した電場だけが小さくなりうる。
:::
::::

::::quiz{answer="-2" tol="1e-6"}
$Z'(0)$ を答えよ。

:::explanation
微分式 $Z'=-2(1+\zeta Z)$ の第 2 項は原点で消える。
:::
::::

::::quiz
正の位相速度で分布の傾きが負、正の波エネルギーの場合に共鳴が与える符号を選べ。

- [x] 減衰
- [ ] 成長
- [ ] 常にゼロ

:::explanation
$\partial_\omega\varepsilon_R>0$、$F'<0$ なら $\varepsilon_I>0$ である。従って $\gamma<0$ となる。
:::
::::

::::quiz{answer="3" tol="1e-6"}
電子の捕捉速度幅について $\Delta\phi$ を九倍にしたときの倍率を答えよ。

:::explanation
分離軌道のエネルギーから速度幅は $\Delta\phi^{1/2}$ に比例し、$\sqrt9=3$ である。
:::
::::
