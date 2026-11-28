---
title: 冷たいプラズマ中の波
description: 線形流体と Maxwell 方程式から、静電波と電磁波の主要な分枝を求める。
notation:
  - symbol: '\omega'
    meaning: 角周波数
  - symbol: '\vb k'
    meaning: 波数ベクトル
  - symbol: 'N'
    meaning: 屈折率
  - symbol: '\Omega_s'
    meaning: 符号付きサイクロトロン角周波数
---

## この章の目標

- 分散行列を作る。
- 電子波とイオン波を区別する。
- 磁場に平行・垂直な波の分枝を求める。

## 線形化と分散関係の求め方

一様平衡について、全摂動を $e^{i(\vb k\cdot\vb x-\omega t)}$ の複素振幅で書き、物理量には実部を用いる。線形性により波数ごとに解ける。背景流がゼロで衝突を無視すれば、:eqref[eq-plasma-species-momentum] と連続式は
$$
-i\omega m_s\vb u_s=q_s(\vb E+\vb u_s\times\vb B_0)-i\vb k\,\gamma_s T_s\,n_{s1}/n_{s0},
\qquad n_{s1}=n_{s0}\vb k\cdot\vb u_s/\omega.
$$
$\gamma_s$ は $\delta p_s/\delta n_s=\gamma_s T_s$ を定める係数である。冷たい模型は $\gamma_s T_s=0$ とする。静電近似は $\vb E=-i\vb k\phi$ として Poisson 方程式を解く。電磁問題では Faraday 則から $\vb B_1=\vb k\times\vb E/\omega$、Ampère–Maxwell 則から
$$
\vb k\times(\vb k\times\vb E)+\frac{\omega^2}{c^2}\vb E
=-i\mu_0\omega\,\vb j.
\label{eq-plasma-wave-equation}
$$
線形速度を $\vb E$ の関数として代入すると同次行列方程式となる。行列式ゼロが非自明解の必要十分条件である。ただし途中で割った $\omega$ や共鳴分母がゼロのときは元の連立式を検討する。

:::example[縦波と横波]
$\vb k=k\hat{\vb x}$ なら二重外積は $(0,-k^2E_y,-k^2E_z)$ である。縦成分の電場は磁場摂動を作らない。
:::

### エネルギー保存と線形応答

冷たい一様プラズマで各種の運動式と $\vb u_s$ の実時間内積を取ると $\partial_t\sum_s m_s n_{s0}u_s^2/2=\vb E\cdot\vb j$ となる。線形 Maxwell 方程式の Poynting 則を足すと
$$
\partial_t\left(\sum_s\frac{m_s n_{s0}u_s^2}{2}+\frac{\varepsilon_0 E^2}{2}+\frac{B_1^2}{2\mu_0}\right)
+\nabla\cdot(\vb E\times\vb B_1/\mu_0)=0.
$$
背景が一様で静止し、外から供給する摂動電流がない場合の式である。周期領域では二次エネルギーの総和が保存され、正の質量密度のもとで指数成長する冷たい線形波はない。異なる平均流を含む場合はこの形の正定値エネルギーをそのまま使えない。

:::example[縦波のエネルギー]
電子だけが動く冷たい縦波の固有周波数は $\omega_{pe}$ である。$u= -ieE/(m_e\omega_{pe})$ とすると平均運動エネルギーは $m_en_0\abs u^2/4=\varepsilon_0\abs E^2/4$ で、平均電場エネルギーと等しい。磁場摂動はゼロである。
:::

### 偏波への射影と初期値の分解

一様・等方な冷たい電子流体でイオンは固定する。任意の波数 $\vb k\ne0$ に対し $\vb k\vb k/k^2$ と $I-\vb k\vb k/k^2$ は直交する射影である。波動行列に作用させると、前者は $(1-\omega_{pe}^2/\omega^2)$、後者は $(1-(\omega_{pe}^2+c^2k^2)/\omega^2)$ の因子を得る。したがって斜めに振れた初期電場を一つの固有波とみなせない。これは空間 Fourier 成分を固定して時間を追う問題である。縦と横の時間初期微分がともにゼロなら、それぞれ固有角周波数の余弦で振動する。縦成分の発散は初期電荷と Gauss の法則で結び、横成分の初期値は電荷に制約されない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $x$ | 波数比 $ck/\omega_{pe}$ | 1 |
| $A$ | 初期電場を規格化した縦成分 | 1 |
| $B$ | 初期電場を規格化した横成分 | 1 |
| $s$ | 時間 $\omega_{pe}t$ | 1 |

横波の規格化角周波数を $\omega_T/\omega_{pe}$ とする。磁場の復元作用が縦波との差になる。

$$
\omega_T/\omega_{pe}=\sqrt{x^{2} + 1}
$$

縦波と横波の位相差を $\Delta\theta$ とする。両者の初期位相をゼロとする。

$$
\Delta\theta=s \left(\omega_T/\omega_{pe} - 1\right)
$$

縦電場を $E_L/E_0$ とする。電荷分離を伴う成分である。

$$
E_L/E_0=A \cos{\left(s \right)}
$$

横電場を $E_T/E_0$ とする。電子慣性と真空の波動項を含む。

$$
E_T/E_0=B \cos{\left(\omega_T/\omega_{pe} s \right)}
$$

電場の大きさを $\abs{\vb E}/E_0$ とする。直交成分を二乗して加える。

$$
\abs{\vb E}/E_0=\sqrt{(E_L/E_0)^2 + (E_T/E_0)^2}
$$

電場エネルギーの縦割合を $f_L$ とする。粒子や磁場のエネルギーは含めない。

$$
f_L=\frac{(E_L/E_0)^2}{(\abs{\vb E}/E_0)^2}
$$

:::example[長波長近くの偏波]
同振幅の二成分が時間とともにどうずれるか計算する。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{10} \\A&=1 \\B&=1 \\s&=\pi
\end{aligned}
$$

これらを順に代入する。

横波の規格化角周波数は

$$
\omega_T/\omega_{pe}=\sqrt{\left(\frac{1}{10}\right)^{2} + 1}\simeq 1.005
$$

である。

縦波と横波の位相差は

$$
\Delta\theta=\left(\pi\right) \left(\left(\frac{\sqrt{101}}{10}\right) - 1\right)\simeq 0.015669
$$

である。

縦電場は

$$
E_L/E_0=\left(1\right) \cos{\left(\left(\pi\right) \right)}=-1
$$

である。

横電場は

$$
E_T/E_0=\left(1\right) \cos{\left(\left(\frac{\sqrt{101}}{10}\right) \left(\pi\right) \right)}\simeq -0.99988
$$

である。

電場の大きさは

$$
\abs{\vb E}/E_0=\sqrt{\left(-1\right)^{2} + \left(\cos{\left(\frac{\sqrt{101} \pi}{10} \right)}\right)^{2}}\simeq 1.4141
$$

である。

電場エネルギーの縦割合は

$$
f_L=\frac{\left(-1\right)^{2}}{\left(\left(\sqrt{\cos^{2}{\left(\frac{\sqrt{101} \pi}{10} \right)} + 1}\right)\right)^{2}}\simeq 0.50006
$$

である。

短い時間では近い角周波数だが、位相差は時間とともに蓄積する。
:::

:::example[片方が節になる時刻]
縦成分が初めてゼロになる時刻を調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=\sqrt{3} \\A&=1 \\B&=1 \\s&=\frac{\pi}{2}
\end{aligned}
$$

これらを順に代入する。

横波の規格化角周波数は

$$
\omega_T/\omega_{pe}=\sqrt{\left(\sqrt{3}\right)^{2} + 1}=2
$$

である。

縦波と横波の位相差は

$$
\Delta\theta=\left(\frac{\pi}{2}\right) \left(\left(2\right) - 1\right)\simeq 1.5708
$$

である。

縦電場は

$$
E_L/E_0=\left(1\right) \cos{\left(\left(\frac{\pi}{2}\right) \right)}=0
$$

である。

横電場は

$$
E_T/E_0=\left(1\right) \cos{\left(\left(2\right) \left(\frac{\pi}{2}\right) \right)}=-1
$$

である。

電場の大きさは

$$
\abs{\vb E}/E_0=\sqrt{\left(-1\right)^{2} + \left(0\right)^{2}}=1
$$

である。

電場エネルギーの縦割合は

$$
f_L=\frac{\left(0\right)^{2}}{\left(1\right)^{2}}=0
$$

である。

横成分だけがその瞬間の電場を担う。
:::

:::example[純粋な電子縦波]
横成分を励起しない初期値を伝播させる。

入力値を次のように置く。

$$
\begin{aligned}
x&=2 \\A&=2 \\B&=0 \\s&=\pi
\end{aligned}
$$

これらを順に代入する。

横波の規格化角周波数は

$$
\omega_T/\omega_{pe}=\sqrt{\left(2\right)^{2} + 1}\simeq 2.2361
$$

である。

縦波と横波の位相差は

$$
\Delta\theta=\left(\pi\right) \left(\left(\sqrt{5}\right) - 1\right)\simeq 3.8832
$$

である。

縦電場は

$$
E_L/E_0=\left(2\right) \cos{\left(\left(\pi\right) \right)}=-2
$$

である。

横電場は

$$
E_T/E_0=\left(0\right) \cos{\left(\left(\pi\right) \left(\sqrt{5}\right) \right)}=0
$$

である。

電場の大きさは

$$
\abs{\vb E}/E_0=\sqrt{\left(-2\right)^{2} + \left(0\right)^{2}}=2
$$

である。

電場エネルギーの縦割合は

$$
f_L=\frac{\left(-2\right)^{2}}{\left(2\right)^{2}}=1
$$

である。

横分枝が存在しても、線形の均一模型では自発的には出現しない。
:::

:::example[純粋な電磁波]
横波の位相を電子縦波の位相と比較する。

入力値を次のように置く。

$$
\begin{aligned}
x&=2 \sqrt{2} \\A&=0 \\B&=1 \\s&=\frac{\pi}{3}
\end{aligned}
$$

これらを順に代入する。

横波の規格化角周波数は

$$
\omega_T/\omega_{pe}=\sqrt{\left(2 \sqrt{2}\right)^{2} + 1}=3
$$

である。

縦波と横波の位相差は

$$
\Delta\theta=\left(\frac{\pi}{3}\right) \left(\left(3\right) - 1\right)\simeq 2.0944
$$

である。

縦電場は

$$
E_L/E_0=\left(0\right) \cos{\left(\left(\frac{\pi}{3}\right) \right)}=0
$$

である。

横電場は

$$
E_T/E_0=\left(1\right) \cos{\left(\left(3\right) \left(\frac{\pi}{3}\right) \right)}=-1
$$

である。

電場の大きさは

$$
\abs{\vb E}/E_0=\sqrt{\left(-1\right)^{2} + \left(0\right)^{2}}=1
$$

である。

電場エネルギーの縦割合は

$$
f_L=\frac{\left(0\right)^{2}}{\left(1\right)^{2}}=0
$$

である。

縦成分に電荷を割り当てない限り、発散のない場が保たれる。
:::


## 電子プラズマ波（Bohm–Gross 波）

非磁化、固定イオンの縦波では、連続式と電子運動量式から $n_{e1}=-n_0e k^2\phi/[m_e(\omega^2-k^2\gamma_e T_e/m_e)]$ となる。この応答を $\varepsilon_0k^2\phi=-en_{e1}$ に代入すると
$$
\omega^2=\omega_{pe}^2+\gamma_e k^2 T_e/m_e.
\label{eq-plasma-bohm-gross}
$$
:term[Bohm–Gross 波]{en="Bohm–Gross wave" yomi="ぼーむぐろすは"} では $\gamma_e=3$ である。波の位相速度が電子熱速度より十分に大きい Maxwell 分布を考える。一次元の無衝突輸送式を線形化する。平衡分布を $f_0(v)=n_0 F(v)$、$\int F\dd v=1$ とすれば、$i(kv-\omega)f_1-e E f_0'/m_e=0$ である。$E=-ik\phi$ を代入して速度積分すると、Poisson 方程式の誘電関数への電子の寄与は $-\omega_{pe}^2 k^{-2}\int F'(v)/(v-\omega/k)\dd v$ となる。ここでは共鳴速度から離れた粒子による実部を漸近展開する。極の因果的な扱いは [運動論](/physics/plasma-physics/kinetic-theory/) で扱う。大きな $\omega/k$ で分母を展開し、部分積分から $\int v F'\dd v=-1$、$\int v^3F'\dd v=-3T_e/m_e$ を使うと、誘電関数は $1-\omega_{pe}^2/\omega^2-3\omega_{pe}^2 k^2T_e/(m_e\omega^4)+\cdots$ となる。高位相速度の根を求めると上の $\gamma_e=3$ と一致する。共鳴粒子による微小な減衰はこの実数の展開では現れない。

:::example[群速度]
$\omega>0$ の枝では $\dv{\omega}{k}=3kv_{te}^2/\omega$ であり、$k\to0$ でゼロになる。冷たい極限では周波数が $k$ に依存しない。
:::

:::example[展開の精度を判定する]
$k\lambda_{De}=0.1$ では Bohm–Gross の温度補正は $\omega^2/\omega_{pe}^2=1.03$、従って周波数の増加は約 $1.49\,\%$ である。対応する位相速度は約 $10.15 v_{te}$ であり、熱速度の大部分から離れている。$k\lambda_{De}=1$ ではこの比は $2$ に下がり、高位相速度の展開と減衰無視の両方を点検する必要がある。
:::

### 二種の可動な冷たい縦波

一様な二種で圧力をゼロにすると、各種の連続式、運動式と Poisson 方程式を合わせて
$$
1-\sum_s\frac{\omega_{ps}^2}{\omega^2}=0
$$
になる。従って非零周波数の電荷振動は一つである。別々の $\omega_{pe}$ と $\omega_{pi}$ を独立な同時固有周波数として足し合わせることはできない。場を駆動する変位は電荷重み付きの和だからである。中性の組み合わせは $\sum_s n_{s0}q_s \xi_s=0$ を満たし、冷たい模型では復元力がない。この零周波数部分は運動式を $\omega^2$ で割って得た誘電式からは消えている。

圧力をともに定数係数 $\gamma_s T_s$ で閉じるなら、$A_s=k^2\gamma_s T_s/m_s$ として
$$
(\omega^2-A_e)(\omega^2-A_i)
=\omega_{pe}^2(\omega^2-A_i)+\omega_{pi}^2(\omega^2-A_e).
$$
これは電子とイオンの変位の二元連立を行列式で直接書いた式でもある。$A_e\ne A_i$ の場合、$\omega^2$ を横軸にした誘電式は二つの極を持つ。各区間での正の微分を使えば、二つの根の一方は両極の間、他方は両極の上にあることが分かる。$A_e=A_i$ なら共通の圧力波が電気的な復元力を受けない固有変位として残る。誘電式の極だけを見る場合には、この縮退を別に調べる。

:::example[等しい圧力応答の縮退]
$A_e=A_i=A$ なら多項式は $(\omega^2-A)(\omega^2-A-\omega_{pe}^2-\omega_{pi}^2)=0$ となる。最初の根では両種の変位を同じにでき、平衡が準中性なら電場が消える。二番目の根では相対変位による復元力が加わる。
:::

### Bohm–Gross 波の位相速度と展開誤差

$v_t^2=T_e/m_e$ とする。短時間の平行圧縮では一次元の無衝突な圧力応答 $\delta p=3T_e\delta n$ により $\omega^2=\omega_{pe}^2+3k^2v_t^2$ となる。流体式の平方根を長波長で展開すれば $1+3(k\lambda_{De})^2/2$ が得られる。平方根の有理化から両者の差は $[3(k\lambda_{De})^2/2]^2/[1+3(k\lambda_{De})^2/2+\sqrt{1+3(k\lambda_{De})^2}]$ である。この差は同じ流体式での展開誤差であり、運動論からの誤差には別途、高次モーメントと共鳴粒子を含める。位相速度が熱速さに近づいたとき、十分に小さい流体展開誤差だけでは減衰の評価ができない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $x$ | $k\lambda_{De}>0$ | 1 |

流体分散式の角周波数比を $\omega/\omega_{pe}$ とする。$\lambda_{De}=v_t/\omega_{pe}$ を入れる。

$$
\omega/\omega_{pe}=\sqrt{3 x^{2} + 1}
$$

二次展開の角周波数比を $\omega_{\rm app}/\omega_{pe}$ とする。$k$ の偶関数なので一次項は消える。

$$
\omega_{\rm app}/\omega_{pe}=\frac{3 x^{2}}{2} + 1
$$

位相速度の比を $v_\phi/v_t$ とする。共鳴速度の尺度でもある。

$$
v_\phi/v_t=\frac{\omega/\omega_{pe}}{x}
$$

群速度の比を $v_g/v_t$ とする。$\omega$ を $k$ で微分する。

$$
v_g/v_t=\frac{3 x}{\omega/\omega_{pe}}
$$

群速度と位相速度の比を $v_g/v_\phi$ とする。電子の復元周波数が有限なので一より小さい。

$$
v_g/v_\phi=\frac{v_g/v_t}{v_\phi/v_t}
$$

展開の相対誤差を $\delta_\omega$ とする。運動論効果はこの比較に含めない。

$$
\delta_\omega=\frac{- \omega/\omega_{pe} + \omega_{\rm app}/\omega_{pe}}{\omega/\omega_{pe}}
$$

:::example[長い波束]
位相速度と群速度を比べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{20}
\end{aligned}
$$

これらを順に代入する。

流体分散式の角周波数比は

$$
\omega/\omega_{pe}=\sqrt{3 \left(\frac{1}{20}\right)^{2} + 1}\simeq 1.0037
$$

である。

二次展開の角周波数比は

$$
\omega_{\rm app}/\omega_{pe}=\frac{3 \left(\frac{1}{20}\right)^{2}}{2} + 1=\frac{803}{800}
$$

である。

位相速度の比は

$$
v_\phi/v_t=\frac{\left(\frac{\sqrt{403}}{20}\right)}{\left(\frac{1}{20}\right)}\simeq 20.075
$$

である。

群速度の比は

$$
v_g/v_t=\frac{3 \left(\frac{1}{20}\right)}{\left(\frac{\sqrt{403}}{20}\right)}\simeq 0.14944
$$

である。

群速度と位相速度の比は

$$
v_g/v_\phi=\frac{\left(\frac{3 \sqrt{403}}{403}\right)}{\left(\sqrt{403}\right)}=\frac{3}{403}
$$

である。

展開の相対誤差は

$$
\delta_\omega=\frac{\left(\frac{803}{800}\right) - \left(\frac{\sqrt{403}}{20}\right)}{\left(\frac{\sqrt{403}}{20}\right)}\simeq 6.9789\times 10^{-6}
$$

である。

位相は速く進んでも狭い波束の包絡が同じ速度で進むわけではない。
:::

:::example[長波長の先頭補正]
周波数の二次補正を流体式と比べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{5}
\end{aligned}
$$

これらを順に代入する。

流体分散式の角周波数比は

$$
\omega/\omega_{pe}=\sqrt{3 \left(\frac{1}{5}\right)^{2} + 1}\simeq 1.0583
$$

である。

二次展開の角周波数比は

$$
\omega_{\rm app}/\omega_{pe}=\frac{3 \left(\frac{1}{5}\right)^{2}}{2} + 1=\frac{53}{50}
$$

である。

位相速度の比は

$$
v_\phi/v_t=\frac{\left(\frac{2 \sqrt{7}}{5}\right)}{\left(\frac{1}{5}\right)}\simeq 5.2915
$$

である。

群速度の比は

$$
v_g/v_t=\frac{3 \left(\frac{1}{5}\right)}{\left(\frac{2 \sqrt{7}}{5}\right)}\simeq 0.56695
$$

である。

群速度と位相速度の比は

$$
v_g/v_\phi=\frac{\left(\frac{3 \sqrt{7}}{14}\right)}{\left(2 \sqrt{7}\right)}=\frac{3}{28}
$$

である。

展開の相対誤差は

$$
\delta_\omega=\frac{- \left(\frac{2 \sqrt{7}}{5}\right) + \left(\frac{53}{50}\right)}{\left(\frac{2 \sqrt{7}}{5}\right)}\simeq 0.0016059
$$

である。

位相速度の熱速さ比も評価して共鳴の尺度を別に把握する。
:::

:::example[Debye 長に迫る波]
流体式と二次展開の差を調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

流体分散式の角周波数比は

$$
\omega/\omega_{pe}=\sqrt{3 \left(\frac{1}{2}\right)^{2} + 1}\simeq 1.3229
$$

である。

二次展開の角周波数比は

$$
\omega_{\rm app}/\omega_{pe}=\frac{3 \left(\frac{1}{2}\right)^{2}}{2} + 1=\frac{11}{8}
$$

である。

位相速度の比は

$$
v_\phi/v_t=\frac{\left(\frac{\sqrt{7}}{2}\right)}{\left(\frac{1}{2}\right)}\simeq 2.6458
$$

である。

群速度の比は

$$
v_g/v_t=\frac{3 \left(\frac{1}{2}\right)}{\left(\frac{\sqrt{7}}{2}\right)}\simeq 1.1339
$$

である。

群速度と位相速度の比は

$$
v_g/v_\phi=\frac{\left(\frac{3 \sqrt{7}}{7}\right)}{\left(\sqrt{7}\right)}=\frac{3}{7}
$$

である。

展開の相対誤差は

$$
\delta_\omega=\frac{\left(\frac{11}{8}\right) - \left(\frac{\sqrt{7}}{2}\right)}{\left(\frac{\sqrt{7}}{2}\right)}\simeq 0.039402
$$

である。

位相速度が高い尾部から中心側へ移り、運動論の検査が重要になる。
:::

:::example[短波長での形式的な極限]
流体式の群速度がどの値に近づくか求める。

入力値を次のように置く。

$$
\begin{aligned}
x&=2
\end{aligned}
$$

これらを順に代入する。

流体分散式の角周波数比は

$$
\omega/\omega_{pe}=\sqrt{3 \left(2\right)^{2} + 1}\simeq 3.6056
$$

である。

二次展開の角周波数比は

$$
\omega_{\rm app}/\omega_{pe}=\frac{3 \left(2\right)^{2}}{2} + 1=7
$$

である。

位相速度の比は

$$
v_\phi/v_t=\frac{\left(\sqrt{13}\right)}{\left(2\right)}\simeq 1.8028
$$

である。

群速度の比は

$$
v_g/v_t=\frac{3 \left(2\right)}{\left(\sqrt{13}\right)}\simeq 1.6641
$$

である。

群速度と位相速度の比は

$$
v_g/v_\phi=\frac{\left(\frac{6 \sqrt{13}}{13}\right)}{\left(\frac{\sqrt{13}}{2}\right)}=\frac{12}{13}
$$

である。

展開の相対誤差は

$$
\delta_\omega=\frac{\left(7\right) - \left(\sqrt{13}\right)}{\left(\sqrt{13}\right)}\simeq 0.94145
$$

である。

これは形式的な流体分枝であり、実在の弱減衰モードを保証しない。
:::


## イオン音波

低周波の縦波では電子慣性を落とし、電子を等温、イオンを冷たいとする。すると $n_{e1}=n_0e\phi/T_e$、$n_{i1}=n_0e k^2\phi/(m_i\omega^2)$ である。Poisson 方程式に代入して
$$
\omega^2=\frac{k^2 c_s^2}{1+k^2\lambda_{De}^2},\qquad c_s=\sqrt{T_e/m_i}
\label{eq-plasma-ion-acoustic}
$$
を得る。これが :term[イオン音波]{en="ion acoustic wave" yomi="いおんおんぱ"} である。有限イオン圧を含めると、イオンの応答の分母が $\omega^2-\gamma_i k^2 T_i/m_i$ になる。そのため分散式は
$$
\omega^2=\gamma_i k^2 T_i/m_i+\frac{k^2 c_s^2}{1+k^2\lambda_{De}^2}
$$
となる。この閉包は $v_{ti}\ll\abs{\omega/k}\ll v_{te}$ で使いやすい。$T_i$ が $T_e$ に近いとイオン共鳴が無視できず、実数の流体根だけでは伝播を判断できない。

:::example[有限 Debye 長の補正]
$k\lambda_{De}=1$ なら位相速度は $c_s/\sqrt2$、群速度は $c_s/(2\sqrt2)$ である。準中性で得る $c_s$ は長波極限の値である。
:::

:::example[有限イオン圧]
等温電子、イオン断熱指数 $\gamma_i=3$、$T_i/T_e=1/30$ の流体模型では、長波長の音速の二乗は $c_s^2(1+3/30)=1.1c_s^2$ である。速度では $\sqrt{1.1}$ 倍となる。運動論でのイオン温度補正の先頭項と一致するが、減衰は流体式に含まれない。
:::

### イオン音波の静電エネルギー

線形化した冷イオン、等温電子の実摂動に対し、電子密度の自由エネルギーは $T_e n_{e1}^2/(2n_0)$ である。イオン運動式から $\partial_t(m_i n_0u_i^2/2)=e n_0 u_i E$。両種の連続式と Gauss 則を使って静電場のエネルギーを加えると、残る仕事率は $e n_0 u_e E$ である。電子の等温つりあいは $e E=-T_e\partial_x n_{e1}/n_0$ だから、これは電子の自由エネルギーの時間微分と発散に分けられる。周期領域での二次総量
$$
W=\int\left(\frac{m_i n_0u_i^2}{2}+\frac{T_e n_{e1}^2}{2n_0}
+\frac{\varepsilon_0E^2}{2}\right)\dd x
$$
は保存される。$k\lambda_{De}\ll1$ では電場エネルギーの比は $k^2\lambda_{De}^2$ と小さく、電子圧力の自由エネルギーが復元力を担う。

:::example[エネルギー分配と短波長]
$E=-ik\phi$、$n_{e1}=n_0 e\phi/T_e$ を入れると、電場エネルギーと電子自由エネルギーの比は $k^2\lambda_{De}^2$ である。$k\lambda_{De}=1/2$ ならこの比は $1/4$ であり、$k\lambda_{De}=2$ なら四倍となる。後者の根を運動論に適用するには、イオンの位相速度が熱速度より大きいという条件を点検する。
:::

### 遮蔽を残したイオン音波の速度

Boltzmann 電子 $\delta n_e=n_0e\phi/T_e$、冷たいイオン $\delta n_i=n_0k^2e\phi/(m_i\omega^2)$ とし、Poisson 方程式に代入する。$k^2\phi=e(\delta n_i-\delta n_e)/\varepsilon_0$ より $\omega^2=k^2c_s^2/(1+k^2\lambda_{De}^2)$ である。密度差は $\delta n_i-\delta n_e=\delta n_e k^2\lambda_{De}^2$ なので、音波の復元に電子圧が効いても電荷分離は小さくできる。群速度は積の微分により $c_s(1+k^2\lambda_{De}^2)^{-3/2}$ となる。高波数では角周波数は $\omega_{pi}$ に漸近するが、有限のイオン温度があればこの近似で分枝を延長できる範囲を再検討する。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $x$ | 電子 Debye 長で規格化した波数 $k\lambda_{De}$ | 1 |
| $a$ | 規格化した電位振幅 $e\phi/T_e$ | 1 |

規格化角周波数を $\omega/\omega_{pi}$ とする。$c_s=\lambda_{De}\omega_{pi}$ を使う。

$$
\omega/\omega_{pi}=\frac{x}{\sqrt{x^{2} + 1}}
$$

位相速度の比を $v_\phi/c_s$ とする。分散式を波数で割る。

$$
v_\phi/c_s=\frac{1}{\sqrt{x^{2} + 1}}
$$

群速度の比を $v_g/c_s$ とする。位相速度の波数依存も含めて微分する。

$$
v_g/c_s=(v_\phi/c_s)^3
$$

イオン密度の摂動を $\delta n_i/n_0$ とする。Poisson 方程式を満たす密度である。

$$
\delta n_i/n_0=a \left(x^{2} + 1\right)
$$

密度差を $(\delta n_i-\delta n_e)/n_0$ とする。電子密度比 $a$ を引く。

$$
(\delta n_i-\delta n_e)/n_0=a x^{2}
$$

イオン摂動に対する電荷分離を $(\delta n_i-\delta n_e)/\delta n_i$ とする。背景に対する比との違いに注意する。

$$
(\delta n_i-\delta n_e)/\delta n_i=\frac{(\delta n_i-\delta n_e)/n_0}{\delta n_i/n_0}
$$

:::example[十分に長い波]
準中性の精度と両速度を比較する。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{20} \\a&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

規格化角周波数は

$$
\omega/\omega_{pi}=\frac{\left(\frac{1}{20}\right)}{\sqrt{\left(\frac{1}{20}\right)^{2} + 1}}\simeq 0.049938
$$

である。

位相速度の比は

$$
v_\phi/c_s=\frac{1}{\sqrt{\left(\frac{1}{20}\right)^{2} + 1}}\simeq 0.99875
$$

である。

群速度の比は

$$
v_g/c_s=\left(\frac{20 \sqrt{401}}{401}\right)^{3}\simeq 0.99626
$$

である。

イオン密度の摂動は

$$
\delta n_i/n_0=\left(\frac{1}{100}\right) \left(\left(\frac{1}{20}\right)^{2} + 1\right)=\frac{401}{40000}
$$

である。

密度差は

$$
(\delta n_i-\delta n_e)/n_0=\left(\frac{1}{100}\right) \left(\frac{1}{20}\right)^{2}=\frac{1}{40000}
$$

である。

イオン摂動に対する電荷分離は

$$
(\delta n_i-\delta n_e)/\delta n_i=\frac{\left(\frac{1}{40000}\right)}{\left(\frac{401}{40000}\right)}=\frac{1}{401}
$$

である。

電子とイオンの摂動はほぼ一致し、波束もほぼ音速で進む。
:::

:::example[中間の波数]
位相速度だけから信号到着を見積もる誤差を調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{2} \\a&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

規格化角周波数は

$$
\omega/\omega_{pi}=\frac{\left(\frac{1}{2}\right)}{\sqrt{\left(\frac{1}{2}\right)^{2} + 1}}\simeq 0.44721
$$

である。

位相速度の比は

$$
v_\phi/c_s=\frac{1}{\sqrt{\left(\frac{1}{2}\right)^{2} + 1}}\simeq 0.89443
$$

である。

群速度の比は

$$
v_g/c_s=\left(\frac{2 \sqrt{5}}{5}\right)^{3}\simeq 0.71554
$$

である。

イオン密度の摂動は

$$
\delta n_i/n_0=\left(\frac{1}{100}\right) \left(\left(\frac{1}{2}\right)^{2} + 1\right)=\frac{1}{80}
$$

である。

密度差は

$$
(\delta n_i-\delta n_e)/n_0=\left(\frac{1}{100}\right) \left(\frac{1}{2}\right)^{2}=\frac{1}{400}
$$

である。

イオン摂動に対する電荷分離は

$$
(\delta n_i-\delta n_e)/\delta n_i=\frac{\left(\frac{1}{400}\right)}{\left(\frac{1}{80}\right)}=\frac{1}{5}
$$

である。

群速度の方が小さく、差は Debye 長による分散から生じる。
:::

:::example[Debye 長の波数]
密度差は摂動の何割か求める。

入力値を次のように置く。

$$
\begin{aligned}
x&=1 \\a&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

規格化角周波数は

$$
\omega/\omega_{pi}=\frac{\left(1\right)}{\sqrt{\left(1\right)^{2} + 1}}\simeq 0.70711
$$

である。

位相速度の比は

$$
v_\phi/c_s=\frac{1}{\sqrt{\left(1\right)^{2} + 1}}\simeq 0.70711
$$

である。

群速度の比は

$$
v_g/c_s=\left(\frac{\sqrt{2}}{2}\right)^{3}\simeq 0.35355
$$

である。

イオン密度の摂動は

$$
\delta n_i/n_0=\left(\frac{1}{100}\right) \left(\left(1\right)^{2} + 1\right)=\frac{1}{50}
$$

である。

密度差は

$$
(\delta n_i-\delta n_e)/n_0=\left(1\right)^{2} \left(\frac{1}{100}\right)=\frac{1}{100}
$$

である。

イオン摂動に対する電荷分離は

$$
(\delta n_i-\delta n_e)/\delta n_i=\frac{\left(\frac{1}{100}\right)}{\left(\frac{1}{50}\right)}=\frac{1}{2}
$$

である。

背景に対して線形でも、摂動同士は準中性で近似しにくい。
:::

:::example[短波長の冷イオン模型]
Poisson 方程式が必要でも振幅が線形かを比べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=4 \\a&=\frac{1}{1000}
\end{aligned}
$$

これらを順に代入する。

規格化角周波数は

$$
\omega/\omega_{pi}=\frac{\left(4\right)}{\sqrt{\left(4\right)^{2} + 1}}\simeq 0.97014
$$

である。

位相速度の比は

$$
v_\phi/c_s=\frac{1}{\sqrt{\left(4\right)^{2} + 1}}\simeq 0.24254
$$

である。

群速度の比は

$$
v_g/c_s=\left(\frac{\sqrt{17}}{17}\right)^{3}\simeq 0.014267
$$

である。

イオン密度の摂動は

$$
\delta n_i/n_0=\left(\frac{1}{1000}\right) \left(\left(4\right)^{2} + 1\right)=\frac{17}{1000}
$$

である。

密度差は

$$
(\delta n_i-\delta n_e)/n_0=\left(4\right)^{2} \left(\frac{1}{1000}\right)=\frac{2}{125}
$$

である。

イオン摂動に対する電荷分離は

$$
(\delta n_i-\delta n_e)/\delta n_i=\frac{\left(\frac{2}{125}\right)}{\left(\frac{17}{1000}\right)}=\frac{16}{17}
$$

である。

イオンの摂動も小さくできるが、電子 Boltzmann 応答と冷イオンの仮定は別に確認する。
:::


## 非磁化プラズマ中の電磁波とカットオフ

冷たい非磁化粒子では $\vb u_s=iq_s\vb E/(m_s\omega)$ である。横偏波に対して :eqref[eq-plasma-wave-equation] に代入すると
$$
\omega^2=c^2k^2+\sum_s\omega_{ps}^2,\qquad
N^2\defeq c^2k^2/\omega^2=1-\sum_s\omega_{ps}^2/\omega^2.
\label{eq-plasma-unmag-em}
$$
$N$ を :term[屈折率]{en="refractive index" yomi="くっせつりつ"} とする。実周波数で $k^2$ が負になる領域では $k$ は純虚数であり、境界から入る波は指数減衰する。$N=0$ は :term[カットオフ]{en="cutoff" yomi="かっとおふ"} である。伝播枝では $v_\phi v_g=c^2$ が微分だけで得られ、エネルギーの移動を表す群速度は $c$ より小さい。

:::example[侵入長]
イオンを無視し $\omega=\omega_{pe}/2$ とする。$1/\abs{k}=2c/(\sqrt3\omega_{pe})$ であり、単純な一様半空間ではこの長さで電場振幅が $e^{-1}$ 倍になる。
:::

### 平面境界での反射

非磁化の一様半空間へ真空から垂直入射する。接線 $E$ と $B$ の連続性から、入射電場を $1$、反射を $r$、透過を $t$ と規格化すると $1+r=t$、$1-r=Nt$ である。従って $r=(1-N)/(1+N)$、$t=2/(1+N)$ となる。実数 $N>0$ なら反射率 $\abs r^2$ と透過率 $N\abs t^2$ の和は $1$ である。$N$ が純虚数なら $\abs r=1$ で、内部の電場は有限であっても平均法線エネルギー流束を運ばない。

:::example[伝播と反射の違い]
$N=1/2$ なら強度の反射率は $1/9$、透過率は $8/9$ である。内部で $k$ が実数であることは、入射エネルギーが全て透過するという意味ではない。
:::

### 密度から反射と遅延を推定する

一様な非磁化プラズマの $\omega^2=\omega_p^2+c^2k^2$ を用いる。$\omega>\omega_p$ なら $N=ck/\omega=\sqrt{1-\omega_p^2/\omega^2}$、真空からの垂直入射の振幅反射係数は $(1-N)/(1+N)$ である。これは電場の接線成分と、磁場の接線成分 $B=\pm NE/c$ を連続にして得られる。群速は $cN$ であり、区間 $L$ の群遅延と真空との差は $(L/c)(1/N-1)$ となる。この値は単一の通過に対する見積もりで、有限層の反射を繰り返すパルスの到着時刻ではない。カットオフより下では $k$ は虚数であり、$N$ を実数と仮定する以下の計算は使わない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $x$ | 密度を臨界密度で割った $\omega_p^2/\omega^2<1$ | 1 |
| $t_0$ | 真空中の区間通過時間 $L/c$ | $\mathrm{s}$ |

屈折率を $N$ とする。実数の正方向の枝を選ぶ。

$$
N=\sqrt{1 - x}
$$

強度反射率を $R$ とする。振幅の二乗である。

$$
R=\frac{\left(1 - N\right)^{2}}{\left(N + 1\right)^{2}}
$$

強度透過率を $T$ とする。透過電場の二乗に媒質の $N$ を掛ける。

$$
T=\frac{4 N}{\left(N + 1\right)^{2}}
$$

位相速度比を $v_\phi/c$ とする。光速を越えても信号の速度ではない。

$$
v_\phi/c=\frac{1}{N}
$$

群速度比を $v_g/c$ とする。分散式を実波数で微分する。

$$
v_g/c=N
$$

真空に対する通過の遅延を $\Delta t$ とする。真空区間の時間を差し引く。

$$
\Delta t=t_{0} \left(-1 + \frac{1}{v_g/c}\right)
$$

:::example[希薄な媒質]
反射と遅延が同じ小パラメータで始まるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{100} \\t_0&=\frac{1}{1000}\,\mathrm{s}
\end{aligned}
$$

これらを順に代入する。

屈折率は

$$
N=\sqrt{1 - \left(\frac{1}{100}\right)}\simeq 0.99499
$$

である。

強度反射率は

$$
R=\frac{\left(1 - \left(\frac{3 \sqrt{11}}{10}\right)\right)^{2}}{\left(\left(\frac{3 \sqrt{11}}{10}\right) + 1\right)^{2}}\simeq 6.3131\times 10^{-6}
$$

である。

強度透過率は

$$
T=\frac{4 \left(\frac{3 \sqrt{11}}{10}\right)}{\left(\left(\frac{3 \sqrt{11}}{10}\right) + 1\right)^{2}}\simeq 0.99999
$$

である。

位相速度比は

$$
v_\phi/c=\frac{1}{\left(\frac{3 \sqrt{11}}{10}\right)}\simeq 1.005
$$

である。

群速度比は

$$
v_g/c=\left(\frac{3 \sqrt{11}}{10}\right)\simeq 0.99499
$$

である。

真空に対する通過の遅延は

$$
\Delta t=\left(\frac{1}{1000}\right) \left(-1 + \frac{1}{\left(\frac{3 \sqrt{11}}{10}\right)}\right)\simeq 5.0378\times 10^{-6}\,\mathrm{s}
$$

である。

反射は密度比の二乗から始まる一方、遅延には密度比の一次がある。
:::

:::example[密度が大きいが伝播する層]
垂直入射で反射電力と通過時間を求める。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{3}{4} \\t_0&=\frac{1}{1000}\,\mathrm{s}
\end{aligned}
$$

これらを順に代入する。

屈折率は

$$
N=\sqrt{1 - \left(\frac{3}{4}\right)}=\frac{1}{2}
$$

である。

強度反射率は

$$
R=\frac{\left(1 - \left(\frac{1}{2}\right)\right)^{2}}{\left(\left(\frac{1}{2}\right) + 1\right)^{2}}=\frac{1}{9}
$$

である。

強度透過率は

$$
T=\frac{4 \left(\frac{1}{2}\right)}{\left(\left(\frac{1}{2}\right) + 1\right)^{2}}=\frac{8}{9}
$$

である。

位相速度比は

$$
v_\phi/c=\frac{1}{\left(\frac{1}{2}\right)}=2
$$

である。

群速度比は

$$
v_g/c=\left(\frac{1}{2}\right)=\frac{1}{2}
$$

である。

真空に対する通過の遅延は

$$
\Delta t=\left(\frac{1}{1000}\right) \left(-1 + \frac{1}{\left(\frac{1}{2}\right)}\right)=\frac{1}{1000}\,\mathrm{s}
$$

である。

群速度は光速の半分であり、位相速度とは逆の傾向を示す。
:::

:::example[カットオフの近傍]
高い反射率と大きな遅延が共存するか調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{99}{100} \\t_0&=\frac{1}{1000}\,\mathrm{s}
\end{aligned}
$$

これらを順に代入する。

屈折率は

$$
N=\sqrt{1 - \left(\frac{99}{100}\right)}=\frac{1}{10}
$$

である。

強度反射率は

$$
R=\frac{\left(1 - \left(\frac{1}{10}\right)\right)^{2}}{\left(\left(\frac{1}{10}\right) + 1\right)^{2}}=\frac{81}{121}
$$

である。

強度透過率は

$$
T=\frac{4 \left(\frac{1}{10}\right)}{\left(\left(\frac{1}{10}\right) + 1\right)^{2}}=\frac{40}{121}
$$

である。

位相速度比は

$$
v_\phi/c=\frac{1}{\left(\frac{1}{10}\right)}=10
$$

である。

群速度比は

$$
v_g/c=\left(\frac{1}{10}\right)=\frac{1}{10}
$$

である。

真空に対する通過の遅延は

$$
\Delta t=\left(\frac{1}{1000}\right) \left(-1 + \frac{1}{\left(\frac{1}{10}\right)}\right)=\frac{9}{1000}\,\mathrm{s}
$$

である。

透過した成分は遅いが入射電力の多くは反射される。
:::

:::example[同じ密度の長区間]
境界での反射と長さに比例する遅延を比べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{3}{4} \\t_0&=\frac{1}{10}\,\mathrm{s}
\end{aligned}
$$

これらを順に代入する。

屈折率は

$$
N=\sqrt{1 - \left(\frac{3}{4}\right)}=\frac{1}{2}
$$

である。

強度反射率は

$$
R=\frac{\left(1 - \left(\frac{1}{2}\right)\right)^{2}}{\left(\left(\frac{1}{2}\right) + 1\right)^{2}}=\frac{1}{9}
$$

である。

強度透過率は

$$
T=\frac{4 \left(\frac{1}{2}\right)}{\left(\left(\frac{1}{2}\right) + 1\right)^{2}}=\frac{8}{9}
$$

である。

位相速度比は

$$
v_\phi/c=\frac{1}{\left(\frac{1}{2}\right)}=2
$$

である。

群速度比は

$$
v_g/c=\left(\frac{1}{2}\right)=\frac{1}{2}
$$

である。

真空に対する通過の遅延は

$$
\Delta t=\left(\frac{1}{10}\right) \left(-1 + \frac{1}{\left(\frac{1}{2}\right)}\right)=\frac{1}{10}\,\mathrm{s}
$$

である。

単一境界の反射率は変わらず、遅延だけが長さに比例する。
:::


## 磁場に垂直な静電波（高域混成波と低域混成波）

$\vb B_0=B_0\hat{\vb z}$、$\vb k=k\hat{\vb x}$ で縦電場を置く。冷たい運動方程式の横二成分を連立すると $u_{sx}=iq_s\omega E_x/[m_s(\omega^2-\Omega_s^2)]$ である。連続式と Poisson 方程式から
$$
S(\omega)=1-\sum_s\frac{\omega_{ps}^2}{\omega^2-\Omega_s^2}=0.
\label{eq-plasma-hybrid}
$$
電子だけ可動なら $\omega_{\rm UH}^2=\omega_{pe}^2+\Omega_e^2$ であり、:term[高域混成波]{en="upper hybrid wave" yomi="こういきこんせいは"} の周波数である。イオンも可動で $\Omega_i^2\ll\omega^2\ll\Omega_e^2$ とすると、$S\simeq1+\omega_{pe}^2/\Omega_e^2-\omega_{pi}^2/(\omega^2-\Omega_i^2)$ となる。したがって :term[低域混成波]{en="lower hybrid wave" yomi="ていいきこんせいは"} は $\omega_{\rm LH}^2\simeq\Omega_i^2+\omega_{pi}^2/(1+\omega_{pe}^2/\Omega_e^2)$ である。準中性の一価プラズマでは $\omega_{pi}^2/\omega_{pe}^2=\Omega_i/\abs{\Omega_e}$ なので、$\omega_{pe}^2\gg\Omega_e^2$ では第 2 項が $\abs{\Omega_i\Omega_e}$ に近づく。純縦電場は有限 $c$ の電磁問題では一般に固有偏波にならず、この式は短波長の静電極限として使う。

:::example[上側の根]
$\omega_{pe}=3\abs{\Omega_e}$ でイオン固定なら $\omega_{\rm UH}=\sqrt{10}\abs{\Omega_e}$ であり、電子プラズマ周波数より上にある。
:::

### 二つの可動種の混成根

$X=\omega^2$、$a_s=\Omega_s^2$、$P_s=\omega_{ps}^2$ と略記する。$a_i<a_e$ に対して $S(X)=1-P_i/(X-a_i)-P_e/(X-a_e)$ である。各極の間で $S'(X)>0$ であり、区間 $(a_i,a_e)$ の両端では $-\infty$ と $+\infty$ に向かう。上側の区間 $(a_e,\infty)$ では $-\infty$ から $1$ まで増加する。従って各区間に一つ根があり、$0\le X<a_i$ では $S>1$ で根がない。分母を払えば二つの根は
$$
(X-a_e)(X-a_i)-P_e(X-a_i)-P_i(X-a_e)=0
$$
で正確に求められる。

:::example[有限のイオン応答]
規格化して $a_i=1,\ a_e=9,\ P_i=1,\ P_e=3$ とすると多項式は $X^2-14X+21$ である。$X=7\pm2\sqrt7$ がそれぞれ二つの区間に入る。実際の電子陽子プラズマでは $a_s$ と $P_s$ は質量比に拘束されるので、この値は分枝の順序を確認するための模型値である。
:::

### 二種の混成周波数を桁落ちなしで求める

垂直静電波で $S=1-\omega_{pe}^2/(\omega^2-\Omega_e^2)-\omega_{pi}^2/(\omega^2-\Omega_i^2)$ とする。$X=\omega^2/\Omega_e^2$ に対して分母を払うと $X^2-AX+C=0$ となる。$r=m_e/m_i$、$p=\omega_{pe}^2/\Omega_e^2$ なら $A=1+r^2+p(1+r)$、$C=r^2+pr^2+pr$ である。小さい根を $(A-\sqrt{A^2-4C})/2$ で求めると、$r$ が小さい場合に桁落ちする。代わりに根の積を用いる。根は極そのものと区別する。$r=1$ の場合には極が一致し、一方の形式根を約分する必要があるので以下は $0<r<1$、$p>0$ に限る。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $r$ | 電子と一価イオンの質量比 | 1 |
| $p$ | 電子プラズマ角周波数と旋回角周波数の二乗比 | 1 |

二根の和を $A$ とする。電子・イオンの両極を払って一次項を整理する。

$$
A=p \left(r + 1\right) + r^{2} + 1
$$

二根の積を $C$ とする。定数項を読む。

$$
C=p r^{2} + p r + r^{2}
$$

上側の二乗根を $X_+$ とする。べきは $X=\omega^2/\Omega_e^2$ に対するものである。

$$
X_+=\frac{A}{2} + \frac{\sqrt{A^{2} - 4 C}}{2}
$$

下側の二乗根を $X_-$ とする。小さい差を避けて計算する。

$$
X_-=\frac{C}{X_{+}}
$$

下側周波数を $\omega_-/\abs{\Omega_e}$ とする。イオン極より上で電子極より下にある。

$$
\omega_-/\abs{\Omega_e}=\sqrt{X_{-}}
$$

電子のみの上側根との比を $X_+/(1+p)$ とする。イオンを固定した高域混成周波数の二乗が $1+p$ である。

$$
X_+/(1+p)=\frac{X_{+}}{p + 1}
$$

:::example[陽子と電子]
二つのスケールの分離を求める。

入力値を次のように置く。

$$
\begin{aligned}
r&=\frac{1}{1836} \\p&=1
\end{aligned}
$$

これらを順に代入する。

二根の和は

$$
A=\left(1\right) \left(\left(\frac{1}{1836}\right) + 1\right) + \left(\frac{1}{1836}\right)^{2} + 1\simeq 2.0005
$$

である。

二根の積は

$$
C=\left(1\right) \left(\frac{1}{1836}\right)^{2} + \left(1\right) \left(\frac{1}{1836}\right) + \left(\frac{1}{1836}\right)^{2}=\frac{919}{1685448}
$$

である。

上側の二乗根は

$$
X_+=\frac{\left(\frac{6743629}{3370896}\right)}{2} + \frac{\sqrt{\left(\frac{6743629}{3370896}\right)^{2} - 4 \left(\frac{919}{1685448}\right)}}{2}\simeq 2.0003
$$

である。

下側の二乗根は

$$
X_-=\frac{\left(\frac{919}{1685448}\right)}{\left(\frac{\sqrt{45451749262249}}{6741792} + \frac{6743629}{6741792}\right)}\simeq 0.00027259
$$

である。

下側周波数は

$$
\omega_-/\abs{\Omega_e}=\sqrt{\left(\frac{919}{1685448 \left(\frac{\sqrt{45451749262249}}{6741792} + \frac{6743629}{6741792}\right)}\right)}\simeq 0.01651
$$

である。

電子のみの上側根との比は

$$
X_+/(1+p)=\frac{\left(\frac{\sqrt{45451749262249}}{6741792} + \frac{6743629}{6741792}\right)}{\left(1\right) + 1}\simeq 1.0001
$$

である。

イオンの周波数が小さくても、下側の根を求めるにはイオンを残す。
:::

:::example[低密度の磁化プラズマ]
下側根がイオン極にどれだけ近いか調べる。

入力値を次のように置く。

$$
\begin{aligned}
r&=\frac{1}{1836} \\p&=\frac{1}{10000}
\end{aligned}
$$

これらを順に代入する。

二根の和は

$$
A=\left(\frac{1}{10000}\right) \left(\left(\frac{1}{1836}\right) + 1\right) + \left(\frac{1}{1836}\right)^{2} + 1\simeq 1.0001
$$

である。

二根の積は

$$
C=\left(\frac{1}{10000}\right) \left(\frac{1}{1836}\right)^{2} + \left(\frac{1}{10000}\right) \left(\frac{1}{1836}\right) + \left(\frac{1}{1836}\right)^{2}\simeq 3.5115\times 10^{-7}
$$

である。

上側の二乗根は

$$
X_+=\frac{\left(\frac{8428085683}{8427240000}\right)}{2} + \frac{\sqrt{- 4 \left(\frac{11837}{33708960000}\right) + \left(\frac{8428085683}{8427240000}\right)^{2}}}{2}\simeq 1.0001
$$

である。

下側の二乗根は

$$
X_-=\frac{\left(\frac{11837}{33708960000}\right)}{\left(\frac{\sqrt{71032528526749696489}}{16854480000} + \frac{8428085683}{16854480000}\right)}\simeq 3.5112\times 10^{-7}
$$

である。

下側周波数は

$$
\omega_-/\abs{\Omega_e}=\sqrt{\left(\frac{11837}{33708960000 \left(\frac{\sqrt{71032528526749696489}}{16854480000} + \frac{8428085683}{16854480000}\right)}\right)}\simeq 0.00059255
$$

である。

電子のみの上側根との比は

$$
X_+/(1+p)=\frac{\left(\frac{\sqrt{71032528526749696489}}{16854480000} + \frac{8428085683}{16854480000}\right)}{\left(\frac{1}{10000}\right) + 1}\simeq 1
$$

である。

$p$ をゼロへ近づけると振動の根は極へ寄るため、密度ゼロで同じモードとはみなせない。
:::

:::example[高密度側]
下側周波数が $\sqrt{r}\abs{\Omega_e}$ に近づくかを調べる。

入力値を次のように置く。

$$
\begin{aligned}
r&=\frac{1}{1836} \\p&=100
\end{aligned}
$$

これらを順に代入する。

二根の和は

$$
A=\left(100\right) \left(\left(\frac{1}{1836}\right) + 1\right) + \left(\frac{1}{1836}\right)^{2} + 1\simeq 101.05
$$

である。

二根の積は

$$
C=\left(100\right) \left(\frac{1}{1836}\right)^{2} + \left(100\right) \left(\frac{1}{1836}\right) + \left(\frac{1}{1836}\right)^{2}=\frac{183701}{3370896}
$$

である。

上側の二乗根は

$$
X_+=\frac{\left(\frac{340644097}{3370896}\right)}{2} + \frac{\sqrt{- 4 \left(\frac{183701}{3370896}\right) + \left(\frac{340644097}{3370896}\right)^{2}}}{2}\simeq 101.05
$$

である。

下側の二乗根は

$$
X_-=\frac{\left(\frac{183701}{3370896}\right)}{\left(\frac{5 \sqrt{4641436954923241}}{6741792} + \frac{340644097}{6741792}\right)}\simeq 0.00053928
$$

である。

下側周波数は

$$
\omega_-/\abs{\Omega_e}=\sqrt{\left(\frac{183701}{3370896 \left(\frac{5 \sqrt{4641436954923241}}{6741792} + \frac{340644097}{6741792}\right)}\right)}\simeq 0.023222
$$

である。

電子のみの上側根との比は

$$
X_+/(1+p)=\frac{\left(\frac{5 \sqrt{4641436954923241}}{6741792} + \frac{340644097}{6741792}\right)}{\left(100\right) + 1}\simeq 1.0005
$$

である。

$p\gg1$ では $X_-\simeq r$ であり、上下で密度依存性が異なる。
:::

:::example[質量差の小さい二種]
質量比を上げて電子だけの根の誤差を評価する。

入力値を次のように置く。

$$
\begin{aligned}
r&=\frac{1}{4} \\p&=2
\end{aligned}
$$

これらを順に代入する。

二根の和は

$$
A=\left(2\right) \left(\left(\frac{1}{4}\right) + 1\right) + \left(\frac{1}{4}\right)^{2} + 1=\frac{57}{16}
$$

である。

二根の積は

$$
C=\left(2\right) \left(\frac{1}{4}\right)^{2} + \left(2\right) \left(\frac{1}{4}\right) + \left(\frac{1}{4}\right)^{2}=\frac{11}{16}
$$

である。

上側の二乗根は

$$
X_+=\frac{\left(\frac{57}{16}\right)}{2} + \frac{\sqrt{- 4 \left(\frac{11}{16}\right) + \left(\frac{57}{16}\right)^{2}}}{2}\simeq 3.3577
$$

である。

下側の二乗根は

$$
X_-=\frac{\left(\frac{11}{16}\right)}{\left(\frac{\sqrt{2545}}{32} + \frac{57}{32}\right)}\simeq 0.20475
$$

である。

下側周波数は

$$
\omega_-/\abs{\Omega_e}=\sqrt{\left(\frac{11}{16 \left(\frac{\sqrt{2545}}{32} + \frac{57}{32}\right)}\right)}\simeq 0.45249
$$

である。

電子のみの上側根との比は

$$
X_+/(1+p)=\frac{\left(\frac{\sqrt{2545}}{32} + \frac{57}{32}\right)}{\left(2\right) + 1}\simeq 1.1192
$$

である。

上側でも可動イオンの寄与が無視しにくくなる。
:::


## 磁場に平行な電磁波（R 波・L 波・ホイッスラー波）

$\vb k\parallel\vb B_0$ なら横成分を円偏波に分けられる。$\vb e_R=(\hat{\vb x}+i\hat{\vb y})/\sqrt2$ は $C\vb v=\vb v\times\hat{\vb z}$ の固有値 $i$、$\vb e_L=\conj{\vb e_R}$ は固有値 $-i$ である。固定点で $\Re(\vb e_R e^{-i\omega t})$ は $+\hat{\vb x}$ から $+\hat{\vb y}$ へ回る。この規約なら
$$
N_R^2=R=1-\sum_s\frac{\omega_{ps}^2}{\omega(\omega+\Omega_s)},\qquad
N_L^2=L=1-\sum_s\frac{\omega_{ps}^2}{\omega(\omega-\Omega_s)} .
\label{eq-plasma-rl}
$$
:term[ホイッスラー波]{en="whistler wave" yomi="ほいっすらーは"} は R 分枝の $\Omega_i\ll\omega\ll\abs{\Omega_e}$ 付近に現れる。イオン項と真空項が小さく $N_R^2\simeq\omega_{pe}^2/(\omega\abs{\Omega_e})$ なら $\omega\simeq\abs{\Omega_e}d_e^2k^2$ である。長い伝播距離では低周波の成分ほど遅れて到着する。

:::example[群速度の倍率]
この二次分散で $k$ を二倍にすると $\omega$ は四倍、$v_g=2\abs{\Omega_e}d_e^2 k$ は二倍である。
:::

### 低周波の両円偏波

準中性なら $\sum_s\omega_{ps}^2/\Omega_s=\sum_s n_s q_s/(\varepsilon_0 B_0)=0$ である。$\omega\ll\abs{\Omega_i},\abs{\Omega_e}$ で :eqref[eq-plasma-rl] を展開すると $1/\omega$ の項が消え、$R,L\simeq1+\sum_s\omega_{ps}^2/\Omega_s^2=1+\rho/(\varepsilon_0B_0^2)$ となる。後者は $1+c^2/v_A^2$ である。$v_A\ll c$ なら両枝とも $\omega/k\simeq v_A$ に近づき、ホイッスラー近似と異なりイオンの応答を落とせない。

:::example[固定イオンとの違い]
電子の項だけでは $\omega\to0$ に $R$ が $1/\omega$ で増大する。一方、可動な準中性イオンが同じ $\omega$ に応答すると逆符号の項を打ち消す。周波数の極限とイオン固定の操作は交換できない。
:::

### Faraday 回転

円偏波の位相速度が異なると、横向きの直線偏波の方向が回転する。両枝が伝播する領域で $E_x$ 偏波を等振幅の R と L の和として $z=0$ から入れる。共通位相を $\bar k z=(k_R+k_L)z/2$、差の半分を $\Delta=(k_R-k_L)z/2$ とすれば
$$
\frac1{\sqrt2}\left(\vb e_R e^{ik_R z}+\vb e_L e^{ik_L z}\right)
=e^{i\bar k z}\left(\hat{\vb x}\cos\Delta-\hat{\vb y}\sin\Delta\right).
$$
この規約で偏波角は $-\Delta$ である。電子項が支配し、$\omega\gg\omega_{pe},\abs{\Omega_e}$ なら $R-L\simeq 2\omega_{pe}^2\Omega_e/\omega^3$ である。平方根を一次まで展開して $k_R-k_L\simeq\omega_{pe}^2\Omega_e/(c\omega^2)$ を得る。従って角度は電子密度と線に沿う磁場の積を距離で積分した量に比例し、周波数の二乗に反比例する。場が連続的に変わる場合は、各円偏波が別々に WKB で進む近似を使う。

:::example[観測周波数を変える]
同じ経路を二倍の周波数で測ると、この高周波極限で偏波角の絶対値は四分の一になる。低周波側が旋回共鳴やカットオフに近い場合は $R,L$ を展開せずに計算する。
:::

### ホイッスラーの曲率と到着順序

イオン応答を無視でき、$\omega\ll\abs{\Omega_e}$ かつ $N^2\gg1$ の R 波では $\omega\simeq\abs{\Omega_e}(kd_e)^2$ である。よって $v_\phi=\abs{\Omega_e}d_e(kd_e)$、$v_g=2v_\phi$ となる。一定区間を各周波数の波束が伝わるとき、時間は $1/\sqrt{\omega}$ に比例する。高周波側が先に着くため、短いインパルスは時刻とともに周波数が下がる信号になる。この関係は同じ経路とほぼ一定の背景を仮定する。実際に $\omega$ を下げ続けるとイオンも動き、低周波で先頭の二乗分散を使えなくなる。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $x$ | 波数 $kd_e$ | 1 |
| $\ell$ | 距離 $L/d_e$ | 1 |
| $r$ | イオンと電子の絶対旋回角周波数の比 | 1 |

近似角周波数を $\omega/\abs{\Omega_e}$ とする。電子の円偏波応答の支配項を残す。

$$
\omega/\abs{\Omega_e}=x^{2}
$$

イオンの旋回に対する比を $\omega/\abs{\Omega_i}$ とする。イオンを固定する順序を確認する。

$$
\omega/\abs{\Omega_i}=\frac{\omega/\abs{\Omega_e}}{r}
$$

位相速度の規格化値を $v_\phi/(\abs{\Omega_e}d_e)$ とする。角周波数を $k$ で割る。

$$
v_\phi/(\abs{\Omega_e}d_e)=x
$$

群速度の規格化値を $v_g/(\abs{\Omega_e}d_e)$ とする。二次の分散を微分する。

$$
v_g/(\abs{\Omega_e}d_e)=2 x
$$

一方向の伝播時間を $t\abs{\Omega_e}$ とする。反射を除いて $L/v_g$ を計算する。

$$
t\abs{\Omega_e}=\frac{\ell}{v_g/(\abs{\Omega_e}d_e)}
$$

四倍の周波数の群遅延比を $t_{\rm next}/t$ とする。四倍でも同じ順序に留まる場合の値である。

$$
t_{\rm next}/t=\frac{1}{2}
$$

:::example[中間の波数]
伝播時間を見積もり、イオンを固定する順序を確かめる。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{10} \\\ell&=100 \\r&=\frac{1}{1836}
\end{aligned}
$$

これらを順に代入する。

近似角周波数は

$$
\omega/\abs{\Omega_e}=\left(\frac{1}{10}\right)^{2}=\frac{1}{100}
$$

である。

イオンの旋回に対する比は

$$
\omega/\abs{\Omega_i}=\frac{\left(\frac{1}{100}\right)}{\left(\frac{1}{1836}\right)}=\frac{459}{25}
$$

である。

位相速度の規格化値は

$$
v_\phi/(\abs{\Omega_e}d_e)=\left(\frac{1}{10}\right)=\frac{1}{10}
$$

である。

群速度の規格化値は

$$
v_g/(\abs{\Omega_e}d_e)=2 \left(\frac{1}{10}\right)=\frac{1}{5}
$$

である。

一方向の伝播時間は

$$
t\abs{\Omega_e}=\frac{\left(100\right)}{\left(\frac{1}{5}\right)}=500
$$

である。

四倍の周波数の群遅延比は

$$
t_{\rm next}/t=\frac{1}{2}
$$

である。

イオン旋回より速いが電子旋回より遅い領域にある。
:::

:::example[二倍の波数]
一つ前の例と到着順序を比べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{5} \\\ell&=100 \\r&=\frac{1}{1836}
\end{aligned}
$$

これらを順に代入する。

近似角周波数は

$$
\omega/\abs{\Omega_e}=\left(\frac{1}{5}\right)^{2}=\frac{1}{25}
$$

である。

イオンの旋回に対する比は

$$
\omega/\abs{\Omega_i}=\frac{\left(\frac{1}{25}\right)}{\left(\frac{1}{1836}\right)}=\frac{1836}{25}
$$

である。

位相速度の規格化値は

$$
v_\phi/(\abs{\Omega_e}d_e)=\left(\frac{1}{5}\right)=\frac{1}{5}
$$

である。

群速度の規格化値は

$$
v_g/(\abs{\Omega_e}d_e)=2 \left(\frac{1}{5}\right)=\frac{2}{5}
$$

である。

一方向の伝播時間は

$$
t\abs{\Omega_e}=\frac{\left(100\right)}{\left(\frac{2}{5}\right)}=250
$$

である。

四倍の周波数の群遅延比は

$$
t_{\rm next}/t=\frac{1}{2}
$$

である。

高い周波数が先に着き、同じ区間の時間は約半分になる。
:::

:::example[イオンに近すぎる波数]
分散式の形式結果と適用条件を照合する。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{100} \\\ell&=100 \\r&=\frac{1}{1836}
\end{aligned}
$$

これらを順に代入する。

近似角周波数は

$$
\omega/\abs{\Omega_e}=\left(\frac{1}{100}\right)^{2}=\frac{1}{10000}
$$

である。

イオンの旋回に対する比は

$$
\omega/\abs{\Omega_i}=\frac{\left(\frac{1}{10000}\right)}{\left(\frac{1}{1836}\right)}=\frac{459}{2500}
$$

である。

位相速度の規格化値は

$$
v_\phi/(\abs{\Omega_e}d_e)=\left(\frac{1}{100}\right)=\frac{1}{100}
$$

である。

群速度の規格化値は

$$
v_g/(\abs{\Omega_e}d_e)=2 \left(\frac{1}{100}\right)=\frac{1}{50}
$$

である。

一方向の伝播時間は

$$
t\abs{\Omega_e}=\frac{\left(100\right)}{\left(\frac{1}{50}\right)}=5000
$$

である。

四倍の周波数の群遅延比は

$$
t_{\rm next}/t=\frac{1}{2}
$$

である。

イオンの比が一より小さいので、この時間を観測予測に使えない。
:::

:::example[電子旋回に接近する波数]
高周波側で近似の前提が保たれるか判定する。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{9}{10} \\\ell&=20 \\r&=\frac{1}{1836}
\end{aligned}
$$

これらを順に代入する。

近似角周波数は

$$
\omega/\abs{\Omega_e}=\left(\frac{9}{10}\right)^{2}=\frac{81}{100}
$$

である。

イオンの旋回に対する比は

$$
\omega/\abs{\Omega_i}=\frac{\left(\frac{81}{100}\right)}{\left(\frac{1}{1836}\right)}=\frac{37179}{25}
$$

である。

位相速度の規格化値は

$$
v_\phi/(\abs{\Omega_e}d_e)=\left(\frac{9}{10}\right)=\frac{9}{10}
$$

である。

群速度の規格化値は

$$
v_g/(\abs{\Omega_e}d_e)=2 \left(\frac{9}{10}\right)=\frac{9}{5}
$$

である。

一方向の伝播時間は

$$
t\abs{\Omega_e}=\frac{\left(20\right)}{\left(\frac{9}{5}\right)}=\frac{100}{9}
$$

である。

四倍の周波数の群遅延比は

$$
t_{\rm next}/t=\frac{1}{2}
$$

である。

$\omega/\abs{\Omega_e}$ が小さくないので R 波の完全な式を解く。
:::


## 磁場に垂直な電磁波（O 波と X 波）

$\vb k=k\hat{\vb x}$ なら $E_z$ の偏波は磁場に平行な電流だけを作る。従って O 波は :eqref[eq-plasma-unmag-em] と同じ分散をもつ。他方、横平面では円偏波を直交成分に戻すと誘電行列は $\begin{pmatrix}S&-iD\\ iD&S\end{pmatrix}$ である。$S=(R+L)/2$、$D=(R-L)/2$ とおく。電磁波の行列はこの $yy$ 成分から $N^2$ を引く。行列式 $S(S-N^2)-D^2=0$ を解いて
$$
N_X^2=RL/S .
\label{eq-plasma-xmode}
$$
:term[異常波]{en="extraordinary wave" yomi="いじょうは"} では $SE_x-iDE_y=0$ であり、通常は縦電場も含む。$S=0$ に近づくとこの模型では大きな $N_X$ が現れ、有限温度や非一様性が重要になる。

:::example[磁場を消す極限]
$\Omega_s=0$ では $R=L$、$D=0$ であり、$N_X^2=S$ になる。O 波との横偏波の分裂は消える。
:::

:::example[偏波と共鳴への接近]
$S\to0$ で $D$ が有限なら、縦行 $SE_x=iDE_y$ は $E_y/E_x=S/(iD)\to0$ を与える。したがって波数が大きくなる X 波の電場は主に $\vb k$ の方向を向く。混成周波数を静電近似から求められる理由である。
:::

### 電子だけが動くときの X 波の通過帯

$\omega>0$ とし、$w_c=\abs{\Omega_e}$、$w_p=\omega_{pe}$ と書く。$R=0$ と $L=0$ の周波数はそれぞれ $\omega_R=(\sqrt{w_c^2+4w_p^2}+w_c)/2$、$\omega_L=(\sqrt{w_c^2+4w_p^2}-w_c)/2$ である。$\omega_R-\omega_L=w_c$ であり、$\omega_R^2-\omega_{\rm UH}^2=w_c(\omega_R-w_c)>0$、$\omega_{\rm UH}^2-\omega_L^2=w_c(w_c+\omega_L)>0$ である。

円偏波の式を通分し、垂直 X 波について
$$
N_X^2=1-\frac{w_p^2(\omega^2-w_p^2)}
{\omega^2(\omega^2-w_c^2-w_p^2)}
=\frac{(\omega^2-\omega_R^2)(\omega^2-\omega_L^2)}
{\omega^2(\omega^2-\omega_{\rm UH}^2)}
$$
を得る。従って $\omega_L<\omega<\omega_{\rm UH}$ と $\omega>\omega_R$ が通過帯である。テンソルの旋回極 $\omega=w_c$ はこの比では約分される。ちょうど極での偏波は、分母で割った運動式によらず元の連立式で扱う。また低い通過帯の周波数がイオンの応答域に入るときには、電子だけを動かす近似を更新する。

### X 波の縦電場と透過窓

電子だけが動く場合に $x=\omega/\abs{\Omega_e}$、$p=\omega_{pe}^2/\Omega_e^2$ とする。Stix の符号規約で $S=1-p/(x^2-1)$、$D=-p/[x(x^2-1)]$ である。垂直伝播を $\hat{\vb x}$ とすれば $S E_x-iD E_y=0$、$(S-N^2)E_y+iDE_x=0$ だから $N_X^2=S-D^2/S$ となる。さらに $E_x=i(D/S)E_y$ である。従って X 波の $N^2$ の零点と、縦電場の比が発散する $S=0$ は別の場所になる。$N_X^2<0$ なら正の実波数は存在せず、境界の近くで減衰する解を選べる。ここでは電子旋回の極と $S=0$ を除く。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $x$ | 電子旋回で規格化した周波数 | 1 |
| $p$ | 電子プラズマ振動数の二乗比 | 1 |

誘電テンソルの対角成分を $S$ とする。種の磁化で周波数依存性が生じる。

$$
S=- \frac{p}{x^{2} - 1} + 1
$$

円偏波の差の半分を $D$ とする。電子の符号付き旋回を使う。

$$
D=- \frac{p}{x \left(x^{2} - 1\right)}
$$

X 波の屈折率の二乗を $N_X^2$ とする。二行の方程式から偏波成分を消去する。

$$
N_X^2=- \frac{D^{2}}{S} + S
$$

O 波の屈折率の二乗を $N_O^2$ とする。電場は平行で旋回を直接駆動しない。

$$
N_O^2=- \frac{p}{x^{2}} + 1
$$

縦成分の振幅比を $\abs{E_x/E_y}$ とする。位相は横成分と四分の一周期ずれる。

$$
\abs{E_x/E_y}=\left|{\frac{D}{S}}\right|
$$

X 波の電場の縦割合を $f_L$ とする。磁場のエネルギーを分母に入れない。

$$
f_L=\frac{\abs{E_x/E_y}^{2}}{\abs{E_x/E_y}^{2} + 1}
$$

:::example[高い周波数での二偏波]
透過できるかと縦成分を求める。

入力値を次のように置く。

$$
\begin{aligned}
x&=3 \\p&=1
\end{aligned}
$$

これらを順に代入する。

誘電テンソルの対角成分は

$$
S=- \frac{\left(1\right)}{\left(3\right)^{2} - 1} + 1=\frac{7}{8}
$$

である。

円偏波の差の半分は

$$
D=- \frac{\left(1\right)}{\left(3\right) \left(\left(3\right)^{2} - 1\right)}=- \frac{1}{24}
$$

である。

X 波の屈折率の二乗は

$$
N_X^2=- \frac{\left(- \frac{1}{24}\right)^{2}}{\left(\frac{7}{8}\right)} + \left(\frac{7}{8}\right)=\frac{55}{63}
$$

である。

O 波の屈折率の二乗は

$$
N_O^2=- \frac{\left(1\right)}{\left(3\right)^{2}} + 1=\frac{8}{9}
$$

である。

縦成分の振幅比は

$$
\abs{E_x/E_y}=\left|{\frac{\left(- \frac{1}{24}\right)}{\left(\frac{7}{8}\right)}}\right|=\frac{1}{21}
$$

である。

X 波の電場の縦割合は

$$
f_L=\frac{\left(\frac{1}{21}\right)^{2}}{\left(\frac{1}{21}\right)^{2} + 1}=\frac{1}{442}
$$

である。

X 波も O 波も伝播でき、高周波で X 波の縦成分は小さい。
:::

:::example[高域混成に近い伝播枝]
$S$ が小さい側で波の向きを調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{7}{5} \\p&=1
\end{aligned}
$$

これらを順に代入する。

誘電テンソルの対角成分は

$$
S=- \frac{\left(1\right)}{\left(\frac{7}{5}\right)^{2} - 1} + 1=- \frac{1}{24}
$$

である。

円偏波の差の半分は

$$
D=- \frac{\left(1\right)}{\left(\frac{7}{5}\right) \left(\left(\frac{7}{5}\right)^{2} - 1\right)}=- \frac{125}{168}
$$

である。

X 波の屈折率の二乗は

$$
N_X^2=- \frac{\left(- \frac{125}{168}\right)^{2}}{\left(- \frac{1}{24}\right)} + \left(- \frac{1}{24}\right)=\frac{649}{49}
$$

である。

O 波の屈折率の二乗は

$$
N_O^2=- \frac{\left(1\right)}{\left(\frac{7}{5}\right)^{2}} + 1=\frac{24}{49}
$$

である。

縦成分の振幅比は

$$
\abs{E_x/E_y}=\left|{\frac{\left(- \frac{125}{168}\right)}{\left(- \frac{1}{24}\right)}}\right|=\frac{125}{7}
$$

である。

X 波の電場の縦割合は

$$
f_L=\frac{\left(\frac{125}{7}\right)^{2}}{\left(\frac{125}{7}\right)^{2} + 1}=\frac{15625}{15674}
$$

である。

X 波の電場は強い縦成分を持ち、O 波との差が明瞭である。
:::

:::example[共鳴の反対側]
共鳴を越えたときにも進行波があるかを調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{3}{2} \\p&=1
\end{aligned}
$$

これらを順に代入する。

誘電テンソルの対角成分は

$$
S=- \frac{\left(1\right)}{\left(\frac{3}{2}\right)^{2} - 1} + 1=\frac{1}{5}
$$

である。

円偏波の差の半分は

$$
D=- \frac{\left(1\right)}{\left(\frac{3}{2}\right) \left(\left(\frac{3}{2}\right)^{2} - 1\right)}=- \frac{8}{15}
$$

である。

X 波の屈折率の二乗は

$$
N_X^2=- \frac{\left(- \frac{8}{15}\right)^{2}}{\left(\frac{1}{5}\right)} + \left(\frac{1}{5}\right)=- \frac{11}{9}
$$

である。

O 波の屈折率の二乗は

$$
N_O^2=- \frac{\left(1\right)}{\left(\frac{3}{2}\right)^{2}} + 1=\frac{5}{9}
$$

である。

縦成分の振幅比は

$$
\abs{E_x/E_y}=\left|{\frac{\left(- \frac{8}{15}\right)}{\left(\frac{1}{5}\right)}}\right|=\frac{8}{3}
$$

である。

X 波の電場の縦割合は

$$
f_L=\frac{\left(\frac{8}{3}\right)^{2}}{\left(\frac{8}{3}\right)^{2} + 1}=\frac{64}{73}
$$

である。

この入力では X 波の二乗屈折率が負、O 波は正である。
:::

:::example[電子プラズマ周波数より低い波]
磁場がある場合の X 波を O 波と比較する。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{3}{4} \\p&=1
\end{aligned}
$$

これらを順に代入する。

誘電テンソルの対角成分は

$$
S=- \frac{\left(1\right)}{\left(\frac{3}{4}\right)^{2} - 1} + 1=\frac{23}{7}
$$

である。

円偏波の差の半分は

$$
D=- \frac{\left(1\right)}{\left(\frac{3}{4}\right) \left(\left(\frac{3}{4}\right)^{2} - 1\right)}=\frac{64}{21}
$$

である。

X 波の屈折率の二乗は

$$
N_X^2=\left(\frac{23}{7}\right) - \frac{\left(\frac{64}{21}\right)^{2}}{\left(\frac{23}{7}\right)}=\frac{95}{207}
$$

である。

O 波の屈折率の二乗は

$$
N_O^2=- \frac{\left(1\right)}{\left(\frac{3}{4}\right)^{2}} + 1=- \frac{7}{9}
$$

である。

縦成分の振幅比は

$$
\abs{E_x/E_y}=\left|{\frac{\left(\frac{64}{21}\right)}{\left(\frac{23}{7}\right)}}\right|=\frac{64}{69}
$$

である。

X 波の電場の縦割合は

$$
f_L=\frac{\left(\frac{64}{69}\right)^{2}}{\left(\frac{64}{69}\right)^{2} + 1}=\frac{4096}{8857}
$$

である。

同じ密度でも X 波に実波数があり、O 波にはない場合がある。
:::


## 演習問題

::::exercise{level=basic}
無磁場の横波で $\omega=2\omega_{pe}$、イオン固定とする。$N^2$ を求めよ。

:::solution
式 :eqref[eq-plasma-unmag-em] から $N^2=1-1/4=3/4$ である。
:::
::::

::::exercise{level=basic}
$\omega_{pe}=\abs{\Omega_e}$ の上部混成周波数を求めよ。

:::solution
電子だけ可動なら二乗を足すので $\omega_{\rm UH}=\sqrt2\omega_{pe}$ である。
:::
::::

::::exercise{level=basic}
電子と共回転する横円偏波を選べ。

:::solution
$\Omega_e<0$ の自由速度は $v_x+iv_y\propto e^{i\abs{\Omega_e}t}$ なので角度が正に増える。固定点の $\vb e_R$ も同じ向きに回る。
:::
::::

::::exercise{level=basic}
$k\lambda_{De}\to0$ のイオン音波の位相速度を求めよ。

:::solution
式 :eqref[eq-plasma-ion-acoustic] で分母が $1$ に近づくため正の枝では $c_s$ に近づく。
:::
::::

::::exercise{level=advanced}
イオン音波の群速度と位相速度の比を導け。

:::solution
$\omega=kc_s(1+k^2\lambda^2)^{-1/2}$ の微分は $c_s(1+k^2\lambda^2)^{-3/2}$ である。$k\ne0$ で位相速度を割ると比は $(1+k^2\lambda^2)^{-1}$ となる。
:::
::::

::::exercise{level=advanced}
Bohm–Gross 波の $\omega/\omega_{pe}$ を $k^2$ の一次まで求めよ。

:::solution
$\omega/\omega_{pe}=(1+3k^2\lambda_{De}^2)^{1/2}$ であり、$k^2$ の一次は $1+(3/2)k^2\lambda_{De}^2$ である。
:::
::::

::::exercise{level=advanced}
一様な非磁化プラズマの横波で、群速度の大きさが $c/2$ となる角周波数を求めよ。

:::solution
$v_g=c^2k/\omega=cN$ なので $N^2=1/4$ である。式 :eqref[eq-plasma-unmag-em] に代入すると $\omega=(2/\sqrt3)\sqrt{\sum_s\omega_{ps}^2}$ である。イオンを固定した場合には和に電子の寄与だけが残る。
:::
::::

::::exercise{level=advanced}
$R=3,L=1$ で X 波の $N^2$ と $E_x/E_y$ を求めよ。

:::solution
$S=2,D=1$ が得られる。$N_X^2=3/2$、縦行の方程式から $E_x/E_y=i/2$ である。
:::
::::

## 確認テスト

::::quiz
O 波の電場は垂直伝播でどの方向を向くか。

- [x] 背景磁場と平行
- [ ] 波数と平行
- [ ] 背景磁場と波数の両方に垂直

:::explanation
$\vb E\parallel\vb B_0$ のとき速度も磁場に平行で磁気力はゼロとなり、通常の横波分散になる。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
$k\lambda_{De}=\sqrt3$ のイオン音波について $\omega/(kc_s)$ を答えよ。

:::explanation
式 :eqref[eq-plasma-ion-acoustic] の正の平方根は $1/\sqrt{1+3}=1/2$ である。
:::
::::

::::quiz
イオンを固定した非磁化プラズマに、真空から $0<\omega<\omega_{pe}$ の横波を垂直入射させると、プラズマ中の波数はどうなるか。

- [x] 純虚数
- [ ] 常に正の実数
- [ ] 常にゼロ

:::explanation
$\omega^2-\omega_{pe}^2<0$ なので $k^2<0$ である。半無限領域では無限遠で減衰する枝を選ぶ。
:::
::::

::::quiz{answer="2.666666667" tol="1e-6"}
$R=4,L=2$ での $N_X^2$ を答えよ。

:::explanation
$S=3$ なので $N_X^2=8/3$ である。
:::
::::
