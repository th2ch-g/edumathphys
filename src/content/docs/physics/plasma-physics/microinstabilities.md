---
title: 速度空間不安定性
description: 速度分布の自由エネルギーが波を成長させる条件を導く。
notation:
  - symbol: 'F'
    meaning: 正規化した平行速度分布
  - symbol: '\gamma'
    meaning: 時間成長率
  - symbol: 'A'
    meaning: 温度異方性 T⊥/T∥−1
---

## この章の目標

- ビーム系の不安定根を求める。
- 誘電関数の実軸上の値から成長を判定する。
- 密度勾配と温度異方性による駆動を比較する。

## 二流体不安定性とビーム–プラズマ不安定性

二つの電子群が反対に流れると、電場が一方からエネルギーを取り出せる。全密度 $n$ を半分ずつに分け、冷たい二群の速度を $\pm U$ とする。固定イオン背景で、各群の連続式と運動式では $\omega$ が $\omega\mp kU$ にずれる。Poisson 方程式に両電荷を足すと
$$
1-\frac{\omega_p^2}{2(\omega-kU)^2}-\frac{\omega_p^2}{2(\omega+kU)^2}=0 .
\label{eq-plasma-two-stream}
$$
$X=\omega^2$、$a=k^2U^2$ とし、分母を払うと $X^2-(2a+\omega_p^2)X+a(a-\omega_p^2)=0$ である。判別式 $\omega_p^4+8a\omega_p^2$ は正である。$0<a<\omega_p^2$ で根の積が負なので、一つの $X$ が負となる。$\omega=i\sqrt{-X}$ は成長根であり、これを :term[二流束不安定性]{en="two-stream instability" yomi="にりゅうそくふあんていせい"} と呼ぶ。

### 希薄ビーム

静止背景の角周波数を $\omega_p$、ビームを $\omega_b$ とする。分散式は $1-\omega_p^2/\omega^2-\omega_b^2/(\omega-kU)^2=0$ である。$kU=\omega_p$ と合わせ $\Delta=\omega-\omega_p$ を小さいとすれば $2\Delta/\omega_p-\omega_b^2/\Delta^2=0$。三乗根のうち上半平面の根の成長率は $(\sqrt3/2)(\omega_p\omega_b^2/2)^{1/3}$ である。ビームの速度幅がこの共鳴幅より広い場合には冷たい極の代わりに速度積分を使う。

:::example[対称ビームの最大値]
負の枝は $X_-=a+\omega_p^2/2-\sqrt{\omega_p^4+8a\omega_p^2}/2$ である。$a$ で微分してゼロとおくと $a=3\omega_p^2/8$、$\gamma_{\max}=\omega_p/\sqrt8$ を得る。
:::

:::example[対称ビームの長波長極限]
$a=k^2U^2\ll\omega_p^2$ なら平方根は $\omega_p^2+4a+\order{a^2/\omega_p^2}$ なので $X_-\simeq-a$ である。従って $\gamma\simeq |kU|$ となる。$k=0$ では分母を払った二次式の零根を元の電場応答に代入できないので、長波長の極限として使う。
:::

### 電子流とイオンの共鳴

同じ密度の電子と一価イオンをともに冷たい流体とし、イオンを静止させる。電子の平均流を $U>0$ とすれば、線形応答は $1-\omega_{pe}^2/(\omega-kU)^2-\omega_{pi}^2/\omega^2=0$ である。この電子・イオン間の流れで起きる成長を :term[Buneman 不安定性]{en="Buneman instability" yomi="ぶねまんふあんていせい"} と呼ぶ。$kU=\omega_{pe}$ の近くでは、イオンの低周波の運動が電子の負の Doppler 周波数の枝と結合する。$kU=\omega_{pe}$、$|\omega|/\omega_{pe}\ll1$ とすれば、電子の応答を展開して
$$
-\frac{2\omega}{\omega_{pe}}-\frac{\omega_{pi}^2}{\omega^2}=0
$$
を得る。$m_e/m_i\ll1$ なら $\omega_{pi}^2/\omega_{pe}^2=m_e/m_i$ である。三根のうち一つは負の実数であり、複素対の成長根は $\gamma=(\sqrt3/2)\omega_{pe}(m_e/2m_i)^{1/3}$ の成長率を持つ。例えば $m_i/m_e=2048$ の模型では $\gamma/\omega_{pe}=\sqrt3/32$ である。電子の速度幅を増すと、冷電子の極は速度積分に置き換わる。電子とイオンの温度を保った電流駆動イオン音波との比較では、この冷流体近似が使えるかを先に調べる。

### 対称な冷ビームの不安定帯を解く

同密度の電子ビーム $\pm U$ と静止正背景を考える。全電子周波数 $\omega_p$ で規格化し $A=(kU/\omega_p)^2$、$X=(\omega/\omega_p)^2$ とする。誘電式 $1-\tfrac12[(\sqrt X-\sqrt A)^{-2}+(\sqrt X+\sqrt A)^{-2}]=0$ から $(X-A)^2-X-A=0$ を得る。従って $X_\pm=(1+2A\pm\sqrt{1+8A})/2$ であり、積 $A(A-1)$ より $0<A<1$ で負根が一つある。成長率は $\sqrt{-X_-}\omega_p$ である。二乗成長率を微分すると $A=3/8$ が最大で、$-X_-=1/8$ になる。$A=0$ のゼロ根は分母の極との合流なので伝播枝から除く。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $A$ | 波数と流速の二乗比 $(kU/\omega_p)^2$ | 1 |

二次式の判別式を $D$ とする。正で二乗周波数が実数となる。

$$
D=8 A + 1
$$

正の分枝を $X_+$ とする。単調な周波数の尺度である。

$$
X_+=A + \frac{\sqrt{D}}{2} + \frac{1}{2}
$$

低い方の二乗周波数を $X_-$ とする。積を使って長波長での桁落ちを防ぐ。

$$
X_-=\frac{A \left(A - 1\right)}{X_{+}}
$$

符号を保持した量 $-X_-$ を求める。正のとき、その平方根が $\gamma/\omega_p$ である。

$$
-X_-=- X_{-}
$$

安定側の振動周波数を $\omega_+/\omega_p$ とする。正の平方根を使う。

$$
\omega_+/\omega_p=\sqrt{X_{+}}
$$

$-8X_-$ は成長帯では最大二乗成長率との比である。安定帯では負になる。

$$
-8X_-=8(-X_-)
$$

:::example[長波長端]
負根の大きさを求める。

入力値を次のように置く。

$$
\begin{aligned}
A&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

二次式の判別式は

$$
D=8 \left(\frac{1}{100}\right) + 1=\frac{27}{25}
$$

である。

正の分枝は

$$
X_+=\left(\frac{1}{100}\right) + \frac{\sqrt{\left(\frac{27}{25}\right)}}{2} + \frac{1}{2}\simeq 1.0296
$$

である。

低い方の二乗周波数は

$$
X_-=\frac{\left(\frac{1}{100}\right) \left(\left(\frac{1}{100}\right) - 1\right)}{\left(\frac{51}{100} + \frac{3 \sqrt{3}}{10}\right)}\simeq -0.0096152
$$

である。

成長率の形式的な二乗は

$$
-X_-=- \left(- \frac{99}{10000 \left(\frac{51}{100} + \frac{3 \sqrt{3}}{10}\right)}\right)\simeq 0.0096152
$$

である。

安定側の振動周波数は

$$
\omega_+/\omega_p=\sqrt{\left(\frac{51}{100} + \frac{3 \sqrt{3}}{10}\right)}\simeq 1.0147
$$

である。

最大二乗成長率との比は

$$
-8X_-=8 \left(\frac{99}{10000 \left(\frac{51}{100} + \frac{3 \sqrt{3}}{10}\right)}\right)\simeq 0.076922
$$

である。

不安定でも $kU$ の小ささで成長率が小さい。
:::

:::example[最速成長]
最大の二乗成長率を計算する。

入力値を次のように置く。

$$
\begin{aligned}
A&=\frac{3}{8}
\end{aligned}
$$

これらを順に代入する。

二次式の判別式は

$$
D=8 \left(\frac{3}{8}\right) + 1=4
$$

である。

正の分枝は

$$
X_+=\frac{\sqrt{\left(4\right)}}{2} + \left(\frac{3}{8}\right) + \frac{1}{2}=\frac{15}{8}
$$

である。

低い方の二乗周波数は

$$
X_-=\frac{\left(\frac{3}{8}\right) \left(\left(\frac{3}{8}\right) - 1\right)}{\left(\frac{15}{8}\right)}=- \frac{1}{8}
$$

である。

成長率の形式的な二乗は

$$
-X_-=- \left(- \frac{1}{8}\right)=\frac{1}{8}
$$

である。

安定側の振動周波数は

$$
\omega_+/\omega_p=\sqrt{\left(\frac{15}{8}\right)}\simeq 1.3693
$$

である。

最大二乗成長率との比は

$$
-8X_-=8 \left(\frac{1}{8}\right)=1
$$

である。

この点で分枝の微分がゼロになる。
:::

:::example[不安定帯の終端]
正根と負根の境界を求める。

入力値を次のように置く。

$$
\begin{aligned}
A&=1
\end{aligned}
$$

これらを順に代入する。

二次式の判別式は

$$
D=8 \left(1\right) + 1=9
$$

である。

正の分枝は

$$
X_+=\left(1\right) + \frac{\sqrt{\left(9\right)}}{2} + \frac{1}{2}=3
$$

である。

低い方の二乗周波数は

$$
X_-=\frac{\left(1\right) \left(\left(1\right) - 1\right)}{\left(3\right)}=0
$$

である。

成長率の形式的な二乗は

$$
-X_-=- \left(0\right)=0
$$

である。

安定側の振動周波数は

$$
\omega_+/\omega_p=\sqrt{\left(3\right)}\simeq 1.7321
$$

である。

最大二乗成長率との比は

$$
-8X_-=8 \left(0\right)=0
$$

である。

成長率はゼロとなるが高い側には振動枝がある。
:::

:::example[短波長の安定帯]
二つの二乗周波数を調べる。

入力値を次のように置く。

$$
\begin{aligned}
A&=2
\end{aligned}
$$

これらを順に代入する。

二次式の判別式は

$$
D=8 \left(2\right) + 1=17
$$

である。

正の分枝は

$$
X_+=\frac{\sqrt{\left(17\right)}}{2} + \left(2\right) + \frac{1}{2}\simeq 4.5616
$$

である。

低い方の二乗周波数は

$$
X_-=\frac{\left(2\right) \left(\left(2\right) - 1\right)}{\left(\frac{\sqrt{17}}{2} + \frac{5}{2}\right)}\simeq 0.43845
$$

である。

成長率の形式的な二乗は

$$
-X_-=- \left(\frac{2}{\frac{\sqrt{17}}{2} + \frac{5}{2}}\right)\simeq -0.43845
$$

である。

安定側の振動周波数は

$$
\omega_+/\omega_p=\sqrt{\left(\frac{\sqrt{17}}{2} + \frac{5}{2}\right)}\simeq 2.1358
$$

である。

最大二乗成長率との比は

$$
-8X_-=8 \left(- \frac{2}{\frac{\sqrt{17}}{2} + \frac{5}{2}}\right)\simeq -3.5076
$$

である。

両方が正なので冷ビームの指数成長はない。
:::


## Penrose の判定条件

無磁化で一様な静電模型では、滑らかな分布の応答を求めるために :eqref[eq-plasma-kinetic-epsilon] を使う。単一種では位相速度を $u=\omega/k$ とする。実軸上で $\Im\varepsilon=-\pi\omega_p^2F'(u)/k^2$ である。$F'(u_0)=0$ かつ十分な減衰を仮定し、部分積分すると
$$
\mathrm{PV}\int\frac{F'(v)}{v-u_0}\dd v
=\int\frac{F(v)-F(u_0)}{(v-u_0)^2}\dd v\equiv P_F(u_0).
$$
境界項は無限遠でゼロで、$u_0$ では $F(v)-F(u_0)=\order{(v-u_0)^2}$ である。従って実軸を横切る位置の $\Re\varepsilon$ は $1-\omega_p^2 P_F(u_0)/k^2$ である。

### 多種の分布をまとめる

運動論で動く種の集合に対し $\omega_*^2=\sum_s\omega_{ps}^2$、$F_{\rm eff}=\sum_s\omega_{ps}^2F_s/\omega_*^2$ とする。式 :eqref[eq-plasma-kinetic-epsilon] の速度積分は $\omega_*^2$ と正規化分布 $F_{\rm eff}$ による一つの積分になる。従って電子とイオンをともに動かす場合は、個々の分布の停留点ではなく $F_{\rm eff}'(v)=0$ となる点を Nyquist 軌跡の実軸交差に使う。固定背景にはこの重みを与えない。以下の証明は $F=F_{\rm eff}$ と読んでも同じである。

### 二峰の場合の証明

二つの極大とその間の極小があり、他に停留点がないとする。さらに両極大で $\Re\varepsilon>0$、いずれの停留点でも $\varepsilon\ne0$ とする。例えば両極大の高さが等しければ、そこで $F(v)-F(u_0)\le0$ なので $\Re\varepsilon\ge1$ である。上半平面で $\varepsilon$ は正則で、大半円では $1$ に近づく。偏角原理により不安定根の個数は、実軸に沿う像の原点周りの回転数である。極小では $F'$ が負から正へ変わるので像が上から下へ実軸を横切る。負の実軸上でこの交差が起これば回転数は $+1$ となる。極大の交差は正の実軸上にあり回転数を変えない。よってこの場合の :term[Penrose の判定条件]{en="Penrose criterion" yomi="ぺんろーずのはんていじょうけん"} は $P_F(u_0)>k^2/\omega_p^2$ である。峰が増える場合は負の実軸を横切る極小と極大の寄与を符号つきで合計する。等号では実軸上の根を別に調べる。

:::example[一峰の分布]
唯一の停留点が全体の最大なら、$F(v)-F(u_0)\le0$ であり $\Re\varepsilon\ge1$ である。負の実軸の交差がなく、この静電模型には上半平面の根がない。
:::

### 水袋分布との照合

速度の区間 $U-d<v<U+d$ で高さ $1/(2d)$ の分布を考える。分布微分は二つの端のデルタの差なので
$$
\int \frac{F'(v)}{v-u}\dd v
=\frac1{2d}\left(\frac1{U-d-u}-\frac1{U+d-u}\right)
=\frac1{(U-u)^2-d^2}.
$$
角のある分布なので滑らかな二峰の証明を直接適用せず、上半平面でこの式を使う。$d\to0$ で冷たい一群の応答に戻る。二つの同じ幅の水袋を重ねると応答の和をとればよく、有限幅の効果は極が各袋の端に移る形で現れる。

:::example[静止水袋の波]
上の式で $U=0$ なら誘電関数の根は $\omega^2=\omega_p^2+k^2d^2$ である。水袋の速度分散は $d^2/3$ なので、Bohm–Gross の $3$ を含む温度補正と一致する。
:::

### 水袋で Penrose 境界の積分を計算する

二つの等密度水袋を $v\in[U-a,U+a]$ と $[-U-a,-U+a]$ に置き、$U>a>0$ とする。全密度で規格化した分布は各区間で $1/(4a)$ である。谷 $v=0$ の分布値はゼロなので、Penrose 型の実数積分 $\int [F(v)-F(0)]/v^2\dd v$ は $1/(U^2-a^2)$ である。速度空間の端では分布が連続でないため、判定を直接滑らかな分布の証明から移す代わりに分散式を求める。各袋の境界デルタ関数を積分すると応答は $-\omega_p^2/[2((\omega\mp kU)^2-k^2a^2)]$ である。$\omega=0$ の零点は $k^2(U^2-a^2)=\omega_p^2$ となる。さらに $\omega=i\gamma$ で二分母の積は正となり、以下のように二乗の根を求められる。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $U_*$ | ビームの中心速度 $kU/\omega_p$ | 1 |
| $A_*$ | 半幅 $ka/\omega_p$ | 1 |

二つの速度尺度の二乗差を $Z$ とする。谷に袋が重ならない入力で正である。

$$
Z=- A_{*}^{2} + U_{*}^{2}
$$

Penrose 積分の波数規格化値を $P$ とする。分布の二つの袋で同じ積分となる。

$$
P=\frac{1}{Z}
$$

二乗周波数の二根の和を $B$ とする。二つの分母の積を展開する。

$$
B=2 A_{*}^{2} + 2 U_{*}^{2} + 1
$$

二根の積を $C$ とする。ゼロ周波数の境界を見る。

$$
C=Z \left(Z - 1\right)
$$

高い二乗根を $X_+$ とする。ここでは正の根である。

$$
X_+=\frac{B}{2} + \frac{\sqrt{B^{2} - 4 C}}{2}
$$

低い二乗根を $X_-$ とする。負なら指数成長がある。

$$
X_-=\frac{C}{X_{+}}
$$

:::example[谷が狭い二峰]
有限幅で積分が十分大きいか調べる。

入力値を次のように置く。

$$
\begin{aligned}
U_*&=1 \\A_*&=\frac{4}{5}
\end{aligned}
$$

これらを順に代入する。

二つの速度尺度の二乗差は

$$
Z=\left(1\right)^{2} - \left(\frac{4}{5}\right)^{2}=\frac{9}{25}
$$

である。

Penrose 積分の波数規格化値は

$$
P=\frac{1}{\left(\frac{9}{25}\right)}=\frac{25}{9}
$$

である。

二乗周波数の二根の和は

$$
B=2 \left(1\right)^{2} + 2 \left(\frac{4}{5}\right)^{2} + 1=\frac{107}{25}
$$

である。

二根の積は

$$
C=\left(\frac{9}{25}\right) \left(\left(\frac{9}{25}\right) - 1\right)=- \frac{144}{625}
$$

である。

高い二乗根は

$$
X_+=\frac{\left(\frac{107}{25}\right)}{2} + \frac{\sqrt{- 4 \left(- \frac{144}{625}\right) + \left(\frac{107}{25}\right)^{2}}}{2}\simeq 4.3332
$$

である。

低い二乗根は

$$
X_-=\frac{\left(- \frac{144}{625}\right)}{\left(\frac{107}{50} + \frac{\sqrt{481}}{10}\right)}\simeq -0.053171
$$

である。

分布の端が谷に近くなると積分が増す。
:::

:::example[ちょうど境界]
谷の積分と低い根を求める。

入力値を次のように置く。

$$
\begin{aligned}
U_*&=\frac{5}{4} \\A_*&=\frac{3}{4}
\end{aligned}
$$

これらを順に代入する。

二つの速度尺度の二乗差は

$$
Z=- \left(\frac{3}{4}\right)^{2} + \left(\frac{5}{4}\right)^{2}=1
$$

である。

Penrose 積分の波数規格化値は

$$
P=\frac{1}{\left(1\right)}=1
$$

である。

二乗周波数の二根の和は

$$
B=2 \left(\frac{3}{4}\right)^{2} + 2 \left(\frac{5}{4}\right)^{2} + 1=\frac{21}{4}
$$

である。

二根の積は

$$
C=\left(1\right) \left(\left(1\right) - 1\right)=0
$$

である。

高い二乗根は

$$
X_+=\frac{\left(\frac{21}{4}\right)}{2} + \frac{\sqrt{- 4 \left(0\right) + \left(\frac{21}{4}\right)^{2}}}{2}=\frac{21}{4}
$$

である。

低い二乗根は

$$
X_-=\frac{\left(0\right)}{\left(\frac{21}{4}\right)}=0
$$

である。

正負の周波数の境界がゼロにある。
:::

:::example[十分離れた袋]
同じ波数で不安定か調べる。

入力値を次のように置く。

$$
\begin{aligned}
U_*&=2 \\A_*&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

二つの速度尺度の二乗差は

$$
Z=\left(2\right)^{2} - \left(\frac{1}{2}\right)^{2}=\frac{15}{4}
$$

である。

Penrose 積分の波数規格化値は

$$
P=\frac{1}{\left(\frac{15}{4}\right)}=\frac{4}{15}
$$

である。

二乗周波数の二根の和は

$$
B=2 \left(2\right)^{2} + 2 \left(\frac{1}{2}\right)^{2} + 1=\frac{19}{2}
$$

である。

二根の積は

$$
C=\left(\frac{15}{4}\right) \left(\left(\frac{15}{4}\right) - 1\right)=\frac{165}{16}
$$

である。

高い二乗根は

$$
X_+=\frac{\left(\frac{19}{2}\right)}{2} + \frac{\sqrt{- 4 \left(\frac{165}{16}\right) + \left(\frac{19}{2}\right)^{2}}}{2}=\frac{33}{4}
$$

である。

低い二乗根は

$$
X_-=\frac{\left(\frac{165}{16}\right)}{\left(\frac{33}{4}\right)}=\frac{5}{4}
$$

である。

谷から離れると積分が小さくなり、負の根は消える。
:::

:::example[冷たい極限に近い幅]
前の冷ビームの結果と比較する。

入力値を次のように置く。

$$
\begin{aligned}
U_*&=\frac{1}{2} \\A_*&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

二つの速度尺度の二乗差は

$$
Z=- \left(\frac{1}{100}\right)^{2} + \left(\frac{1}{2}\right)^{2}=\frac{2499}{10000}
$$

である。

Penrose 積分の波数規格化値は

$$
P=\frac{1}{\left(\frac{2499}{10000}\right)}=\frac{10000}{2499}
$$

である。

二乗周波数の二根の和は

$$
B=2 \left(\frac{1}{100}\right)^{2} + 2 \left(\frac{1}{2}\right)^{2} + 1=\frac{7501}{5000}
$$

である。

二根の積は

$$
C=\left(\frac{2499}{10000}\right) \left(\left(\frac{2499}{10000}\right) - 1\right)\simeq -0.18745
$$

である。

高い二乗根は

$$
X_+=\frac{\left(\frac{7501}{5000}\right)}{2} + \frac{\sqrt{- 4 \left(- \frac{18744999}{100000000}\right) + \left(\frac{7501}{5000}\right)^{2}}}{2}\simeq 1.6162
$$

である。

低い二乗根は

$$
X_-=\frac{\left(- \frac{18744999}{100000000}\right)}{\left(\frac{7501}{10000} + \frac{\sqrt{7501}}{100}\right)}\simeq -0.11598
$$

である。

幅の影響が小さく冷たい二流に近い。
:::


### 滑らかな二峰分布の複素根

二つの等密度ビームに対し、単位密度あたりの分布を $F_\pm(v)=(\Delta/\pi)[(v\mp U)^2+\Delta^2]^{-1}$ とする。$\Delta>0$ は半値半幅である。$k>0$ とし、$\Im u>0$ で部分積分して $\int F_\pm'(v)/(v-u)\dd v=\int F_\pm(v)/(v-u)^2\dd v$ を得る。下半平面に閉じると半円の寄与は消え、時計回りなので $-2\pi i$ と下側の極の残差を掛ける。その結果は $(u\mp U+i\Delta)^{-2}$ である。従って $\tilde\omega=\omega+ik\Delta$ を用いると冷ビームと同じ式 :eqref[eq-plasma-two-stream] になる。ただしこの変換は上半平面に根があるかの判定に使う。$A=(kU/\omega_p)^2<1$ の冷たい負の二乗根 $X_-$ に対して成長率は $\gamma/\omega_p=\sqrt{-X_-}-k\Delta/\omega_p$ である。二峰の広がりは成長を弱め、冷たいビームでの不安定帯全体を保つとは限らない。

この分布は速度の尾部が $v^{-2}$ なので二次モーメントが発散する。ここでは滑らかさを保ったまま誘電関数の積分を解析的に解くために使っている。$\Delta$ を Maxwell 分布の温度へ読み替えない。分布の中心 $v=0$ の二階微分は $2\Delta(3U^2-\Delta^2)/[\pi(U^2+\Delta^2)^3]$ であり、$\Delta<\sqrt3U$ なら谷があるが、その条件だけでは与えた波数の不安定性を保証しない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $A$ | 波数比の二乗 $(kU/\omega_p)^2$。$0<A<1$ とする。 | 1 |
| $d$ | ビーム幅と流速の比 $\Delta/U$ | 1 |

高い方の冷たい二乗根を $X_+$ とする。冷たい分散式を用いる。

$$
X_+=A + \frac{\sqrt{8 A + 1}}{2} + \frac{1}{2}
$$

低い方の冷たい二乗根を $X_-$ とする。負なので上半平面に候補がある。

$$
X_-=\frac{A \left(A - 1\right)}{X_{+}}
$$

冷たい場合の成長率を $G_0$ とする。単位は全電子の $\omega_p$ である。

$$
G_0=\sqrt{- X_{-}}
$$

速度幅の補正を $G_\Delta$ とする。波数と速度幅を掛けて周波数にする。

$$
G_\Delta=\sqrt{A} d
$$

有限幅を含む符号付き成長率の候補を $G$ とする。正なら上半平面に零点が残る。

$$
G=G_{0} - G_{\Delta}
$$

その波数でのしきい幅を $d_c$ とする。冷たい成長率と速度幅が等しい。

$$
d_c=\frac{G_{0}}{\sqrt{A}}
$$

:::example[狭いビーム]
冷たい波数の最大位置で幅の効果を求める。

入力値を次のように置く。

$$
\begin{aligned}
A&=\frac{3}{8} \\d&=\frac{1}{10}
\end{aligned}
$$

これらを順に代入する。

高い方の冷たい二乗根は

$$
X_+=\left(\frac{3}{8}\right) + \frac{\sqrt{8 \left(\frac{3}{8}\right) + 1}}{2} + \frac{1}{2}=\frac{15}{8}
$$

である。

低い方の冷たい二乗根は

$$
X_-=\frac{\left(\frac{3}{8}\right) \left(\left(\frac{3}{8}\right) - 1\right)}{\left(\frac{15}{8}\right)}=- \frac{1}{8}
$$

である。

冷たい場合の成長率は

$$
G_0=\sqrt{- \left(- \frac{1}{8}\right)}\simeq 0.35355
$$

である。

速度幅の補正は

$$
G_\Delta=\left(\frac{1}{10}\right) \sqrt{\left(\frac{3}{8}\right)}\simeq 0.061237
$$

である。

有限幅を含む符号付き成長率の候補は

$$
G=\left(\frac{\sqrt{2}}{4}\right) - \left(\frac{\sqrt{6}}{40}\right)\simeq 0.29232
$$

である。

その波数でのしきい幅は

$$
d_c=\frac{\left(\frac{\sqrt{2}}{4}\right)}{\sqrt{\left(\frac{3}{8}\right)}}\simeq 0.57735
$$

である。

幅は成長を抑えるが、まだ正の根が残る。
:::

:::example[しきいに近い幅]
谷があるまま成長がどの程度減るか計算する。

入力値を次のように置く。

$$
\begin{aligned}
A&=\frac{3}{8} \\d&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

高い方の冷たい二乗根は

$$
X_+=\left(\frac{3}{8}\right) + \frac{\sqrt{8 \left(\frac{3}{8}\right) + 1}}{2} + \frac{1}{2}=\frac{15}{8}
$$

である。

低い方の冷たい二乗根は

$$
X_-=\frac{\left(\frac{3}{8}\right) \left(\left(\frac{3}{8}\right) - 1\right)}{\left(\frac{15}{8}\right)}=- \frac{1}{8}
$$

である。

冷たい場合の成長率は

$$
G_0=\sqrt{- \left(- \frac{1}{8}\right)}\simeq 0.35355
$$

である。

速度幅の補正は

$$
G_\Delta=\left(\frac{1}{2}\right) \sqrt{\left(\frac{3}{8}\right)}\simeq 0.30619
$$

である。

有限幅を含む符号付き成長率の候補は

$$
G=\left(\frac{\sqrt{2}}{4}\right) - \left(\frac{\sqrt{6}}{8}\right)\simeq 0.047367
$$

である。

その波数でのしきい幅は

$$
d_c=\frac{\left(\frac{\sqrt{2}}{4}\right)}{\sqrt{\left(\frac{3}{8}\right)}}\simeq 0.57735
$$

である。

小さな正根が残る。
:::

:::example[谷があっても成長しない波数]
滑らかな二峰という条件だけで不安定かを調べる。

入力値を次のように置く。

$$
\begin{aligned}
A&=\frac{3}{8} \\d&=\frac{4}{5}
\end{aligned}
$$

これらを順に代入する。

高い方の冷たい二乗根は

$$
X_+=\left(\frac{3}{8}\right) + \frac{\sqrt{8 \left(\frac{3}{8}\right) + 1}}{2} + \frac{1}{2}=\frac{15}{8}
$$

である。

低い方の冷たい二乗根は

$$
X_-=\frac{\left(\frac{3}{8}\right) \left(\left(\frac{3}{8}\right) - 1\right)}{\left(\frac{15}{8}\right)}=- \frac{1}{8}
$$

である。

冷たい場合の成長率は

$$
G_0=\sqrt{- \left(- \frac{1}{8}\right)}\simeq 0.35355
$$

である。

速度幅の補正は

$$
G_\Delta=\sqrt{\left(\frac{3}{8}\right)} \left(\frac{4}{5}\right)\simeq 0.4899
$$

である。

有限幅を含む符号付き成長率の候補は

$$
G=\left(\frac{\sqrt{2}}{4}\right) - \left(\frac{\sqrt{6}}{5}\right)\simeq -0.13634
$$

である。

その波数でのしきい幅は

$$
d_c=\frac{\left(\frac{\sqrt{2}}{4}\right)}{\sqrt{\left(\frac{3}{8}\right)}}\simeq 0.57735
$$

である。

形式根が上半平面にない。幅がしきいを越えている。
:::

:::example[長波長での有限幅]
同じ分布で波数を減らした場合を調べる。

入力値を次のように置く。

$$
\begin{aligned}
A&=\frac{1}{100} \\d&=\frac{4}{5}
\end{aligned}
$$

これらを順に代入する。

高い方の冷たい二乗根は

$$
X_+=\left(\frac{1}{100}\right) + \frac{\sqrt{8 \left(\frac{1}{100}\right) + 1}}{2} + \frac{1}{2}\simeq 1.0296
$$

である。

低い方の冷たい二乗根は

$$
X_-=\frac{\left(\frac{1}{100}\right) \left(\left(\frac{1}{100}\right) - 1\right)}{\left(\frac{51}{100} + \frac{3 \sqrt{3}}{10}\right)}\simeq -0.0096152
$$

である。

冷たい場合の成長率は

$$
G_0=\sqrt{- \left(- \frac{99}{10000 \left(\frac{51}{100} + \frac{3 \sqrt{3}}{10}\right)}\right)}\simeq 0.098057
$$

である。

速度幅の補正は

$$
G_\Delta=\sqrt{\left(\frac{1}{100}\right)} \left(\frac{4}{5}\right)=\frac{2}{25}
$$

である。

有限幅を含む符号付き成長率の候補は

$$
G=- \left(\frac{2}{25}\right) + \left(\frac{3 \sqrt{11}}{100 \sqrt{\frac{51}{100} + \frac{3 \sqrt{3}}{10}}}\right)\simeq 0.018057
$$

である。

その波数でのしきい幅は

$$
d_c=\frac{\left(\frac{3 \sqrt{11}}{100 \sqrt{\frac{51}{100} + \frac{3 \sqrt{3}}{10}}}\right)}{\sqrt{\left(\frac{1}{100}\right)}}\simeq 0.98057
$$

である。

幅の Doppler 項も小さくなるので、前例で安定でもこの波数では成長が残る。
:::


## イオン音波不安定性

イオン音波では静止イオンが波からエネルギーを受け取りやすい。電子の平均速度 $U_e$ が正なら、電子の分布の傾きは $v_\phi-U_e$ の符号に依存する。電子をずれた Maxwell 分布、イオンを静止 Maxwell 分布とし、$\abs{v_\phi-U_e}\ll v_{te}$、$v_\phi\gg v_{ti}$ とする。式 :eqref[eq-plasma-landau-rate] で $\varepsilon_I$ の共通正因子 $\sqrt{\pi/2}/(k^2\lambda_{De}^2)$ を除くと
$$
\frac{v_\phi-U_e}{v_{te}}+
\frac{T_e}{T_i}\frac{v_\phi}{v_{ti}}e^{-v_\phi^2/(2v_{ti}^2)}
$$
となる。正の枝で $\partial_\omega\varepsilon_R>0$ なら、これが負のとき成長する。電子の流れが位相速度を越えただけでは、有限イオン温度での減衰を必ず上回るとは限らない。

:::example[冷たいイオン近似]
$T_i/T_e\to0$ では指数因子がべき因子に勝つ。従ってしきい値は $U_e=v_\phi$ に近づく。
:::

:::example[流れを増やすときの符号]
$U_e=2v_\phi$、$v_\phi/v_{te}=0.02$ なら電子の寄与は $-0.02$ である。イオン項が $0.005$ なら合計が負となり成長するが、$0.03$ なら合計が正となり減衰する。この比較は同じ実位相速度と弱い虚部を仮定している。
:::

### イオン音波に流れが与える逆符号の寄与

電子の平均流を $U_e$、イオンの平均流をゼロとする。位相速度 $v_\phi>0$ で $\zeta_e=(v_\phi-U_e)/(\sqrt2v_{te})$、$\zeta_i=v_\phi/(\sqrt2v_{ti})$ である。Maxwell 応答の虚部を共通の $\sqrt\pi/(k\lambda_{De})^2$ で割ると $G=\zeta_e e^{-\zeta_e^2}+\tau^{-1}\zeta_i e^{-\zeta_i^2}$ である。正エネルギーで弱減衰の分枝なら $G<0$ が電子流による成長の指標になる。$U_e>v_\phi$ は電子寄与を反転させるが、全体に対しては十分条件ではない。また大きな $U_e$ で実部の零点が移るため、以下では位相速度を指定して虚部だけを比較する。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $x$ | 電子から見た共鳴引数 $\zeta_e$ | 1 |
| $y$ | イオンから見た共鳴引数 $\zeta_i$ | 1 |
| $\tau$ | 温度比 $T_i/T_e$ | 1 |

電子の符号付き寄与を $G_e$ とする。流れが位相速度を超えると負になる。

$$
G_e=x e^{- x^{2}}
$$

イオンの符号付き寄与を $G_i$ とする。静止イオンはここでは吸収する。

$$
G_i=\frac{y e^{- y^{2}}}{\tau}
$$

全体の虚部因子を $G$ とする。正なら弱減衰側である。

$$
G=G_{e} + G_{i}
$$

電子とイオンの大きさの比を $\abs{G_e}/G_i$ とする。電子が負の場合は一と比較する。

$$
\abs{G_e}/G_i=\frac{\left|{G_{e}}\right|}{G_{i}}
$$

イオン共鳴での分布比を $F_i/F_{i0}$ とする。温度比とは別の抑制要因である。

$$
F_i/F_{i0}=e^{- y^{2}}
$$

最大電子寄与に対する比を $\sqrt{2e}\abs{G_e}$ とする。$|x|=1/\sqrt2$ で最大になる。

$$
\sqrt{2e}\abs{G_e}=\sqrt{2} e^{\frac{1}{2}} \left|{G_{e}}\right|
$$

:::example[電子流が遅い場合]
両種の符号を比較する。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{20} \\y&=4 \\\tau&=\frac{1}{30}
\end{aligned}
$$

これらを順に代入する。

電子の符号付き寄与は

$$
G_e=\left(\frac{1}{20}\right) e^{- \left(\frac{1}{20}\right)^{2}}\simeq 0.049875
$$

である。

イオンの符号付き寄与は

$$
G_i=\frac{\left(4\right) e^{- \left(4\right)^{2}}}{\left(\frac{1}{30}\right)}\simeq 1.3504\times 10^{-5}
$$

である。

全体の虚部因子は

$$
G=\left(\frac{120}{e^{16}}\right) + \left(\frac{1}{20 e^{\frac{1}{400}}}\right)\simeq 0.049889
$$

である。

電子とイオンの大きさの比は

$$
\abs{G_e}/G_i=\frac{\left|{\left(\frac{1}{20 e^{\frac{1}{400}}}\right)}\right|}{\left(\frac{120}{e^{16}}\right)}\simeq 3693.3
$$

である。

イオン共鳴での分布比は

$$
F_i/F_{i0}=e^{- \left(4\right)^{2}}\simeq 1.1254\times 10^{-7}
$$

である。

最大電子寄与に対する比は

$$
\sqrt{2e}\abs{G_e}=\sqrt{2} e^{\frac{1}{2}} \left|{\left(\frac{1}{20 e^{\frac{1}{400}}}\right)}\right|\simeq 0.11629
$$

である。

電子もイオンも吸収する。
:::

:::example[電子流による供給]
電子を位相速度より速くし全寄与を調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=- \frac{1}{20} \\y&=4 \\\tau&=\frac{1}{30}
\end{aligned}
$$

これらを順に代入する。

電子の符号付き寄与は

$$
G_e=\left(- \frac{1}{20}\right) e^{- \left(- \frac{1}{20}\right)^{2}}\simeq -0.049875
$$

である。

イオンの符号付き寄与は

$$
G_i=\frac{\left(4\right) e^{- \left(4\right)^{2}}}{\left(\frac{1}{30}\right)}\simeq 1.3504\times 10^{-5}
$$

である。

全体の虚部因子は

$$
G=\left(- \frac{1}{20 e^{\frac{1}{400}}}\right) + \left(\frac{120}{e^{16}}\right)\simeq -0.049862
$$

である。

電子とイオンの大きさの比は

$$
\abs{G_e}/G_i=\frac{\left|{\left(- \frac{1}{20 e^{\frac{1}{400}}}\right)}\right|}{\left(\frac{120}{e^{16}}\right)}\simeq 3693.3
$$

である。

イオン共鳴での分布比は

$$
F_i/F_{i0}=e^{- \left(4\right)^{2}}\simeq 1.1254\times 10^{-7}
$$

である。

最大電子寄与に対する比は

$$
\sqrt{2e}\abs{G_e}=\sqrt{2} e^{\frac{1}{2}} \left|{\left(- \frac{1}{20 e^{\frac{1}{400}}}\right)}\right|\simeq 0.11629
$$

である。

冷たいイオンの寄与が小さく電子の供給が勝つ。
:::

:::example[温かいイオンによる抑止]
電子の符号を保ってイオンの吸収を変える。

入力値を次のように置く。

$$
\begin{aligned}
x&=- \frac{1}{20} \\y&=2 \\\tau&=\frac{1}{5}
\end{aligned}
$$

これらを順に代入する。

電子の符号付き寄与は

$$
G_e=\left(- \frac{1}{20}\right) e^{- \left(- \frac{1}{20}\right)^{2}}\simeq -0.049875
$$

である。

イオンの符号付き寄与は

$$
G_i=\frac{\left(2\right) e^{- \left(2\right)^{2}}}{\left(\frac{1}{5}\right)}\simeq 0.18316
$$

である。

全体の虚部因子は

$$
G=\left(- \frac{1}{20 e^{\frac{1}{400}}}\right) + \left(\frac{10}{e^{4}}\right)\simeq 0.13328
$$

である。

電子とイオンの大きさの比は

$$
\abs{G_e}/G_i=\frac{\left|{\left(- \frac{1}{20 e^{\frac{1}{400}}}\right)}\right|}{\left(\frac{10}{e^{4}}\right)}\simeq 0.27231
$$

である。

イオン共鳴での分布比は

$$
F_i/F_{i0}=e^{- \left(2\right)^{2}}\simeq 0.018316
$$

である。

最大電子寄与に対する比は

$$
\sqrt{2e}\abs{G_e}=\sqrt{2} e^{\frac{1}{2}} \left|{\left(- \frac{1}{20 e^{\frac{1}{400}}}\right)}\right|\simeq 0.11629
$$

である。

流れが位相速度を超えていても吸収が勝る。
:::

:::example[電子の中心から遠い共鳴]
速い流れで寄与が単調に増すか調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=-3 \\y&=4 \\\tau&=\frac{1}{30}
\end{aligned}
$$

これらを順に代入する。

電子の符号付き寄与は

$$
G_e=\left(-3\right) e^{- \left(-3\right)^{2}}\simeq -0.00037023
$$

である。

イオンの符号付き寄与は

$$
G_i=\frac{\left(4\right) e^{- \left(4\right)^{2}}}{\left(\frac{1}{30}\right)}\simeq 1.3504\times 10^{-5}
$$

である。

全体の虚部因子は

$$
G=\left(- \frac{3}{e^{9}}\right) + \left(\frac{120}{e^{16}}\right)\simeq -0.00035673
$$

である。

電子とイオンの大きさの比は

$$
\abs{G_e}/G_i=\frac{\left|{\left(- \frac{3}{e^{9}}\right)}\right|}{\left(\frac{120}{e^{16}}\right)}\simeq 27.416
$$

である。

イオン共鳴での分布比は

$$
F_i/F_{i0}=e^{- \left(4\right)^{2}}\simeq 1.1254\times 10^{-7}
$$

である。

最大電子寄与に対する比は

$$
\sqrt{2e}\abs{G_e}=\sqrt{2} e^{\frac{1}{2}} \left|{\left(- \frac{3}{e^{9}}\right)}\right|\simeq 0.00086324
$$

である。

電子の共鳴位置が尾部へ移るので寄与は再び弱くなる。
:::


## ドリフト波と普遍的不安定性

一定の磁場 $B\hat{\vb z}$ と、$L_n^{-1}=-\partial_x\ln n_0>0$ の密度勾配を考える。局所近似で背景の係数を固定し、冷たいイオン、断熱電子 $\delta n_e/n_0=e\phi/T_e$、無視できるイオン平行運動を仮定する。各種の仮定を同時に使うには $|k_\parallel|c_s\ll|\omega|\ll |k_\parallel|v_{te}$ が必要である。イオンの連続式に $\vb v_E=(-\partial_y\phi,\partial_x\phi,0)/B$ と $\vb v_{\rm pol}=-m_i\partial_t\nabla_\perp\phi/(eB^2)$ を入れると
$$
-i\omega (1+k_\perp^2\rho_s^2)\frac{e\phi}{T_e}
+i\frac{k_y}{BL_n}\phi=0,\qquad
\rho_s^2=\frac{T_e m_i}{e^2B^2}.
$$
よって :term[ドリフト波]{en="drift wave" yomi="どりふとは"} は $\omega=\omega_{*e}/(1+k_\perp^2\rho_s^2)$、$\omega_{*e}=k_y T_e/(eBL_n)$ である。この閉包では根は実数である。

### 電子の位相差

電子平行の速度分布が Maxwell 型で、密度勾配を保ち、旋回半径を無視する。線形ドリフト運動論では
$$
(\partial_t+v_\parallel\partial_z)f_{e1}
=-v_{E,x}\partial_x f_{e0}-(q_e E_\parallel/m_e)\partial_{v_\parallel}f_{e0}.
$$
これを Fourier 化し、右辺の両項を合わせると $f_{e1}=(e\phi/T_e)f_{e0}[1-(\omega-\omega_{*e})/(\omega-k_\parallel v_\parallel)]$ となる。$k_\parallel>0$ として速度積分をとると $\delta n_e/n_0=(e\phi/T_e)[1+(\omega-\omega_{*e})Z(\zeta_e)/(\sqrt2 k_\parallel v_{te})]$。イオン密度は $(e\phi/T_e)[\omega_{*e}/\omega-k_\perp^2\rho_s^2]$ である。小さい電子共鳴補正を入れると、分散式の虚部は $\sqrt\pi(\omega_r-\omega_{*e})/(\sqrt2 k_\parallel v_{te})$ である。$|\omega_r|\ll k_\parallel v_{te}$ かつ $0<\omega_r<\omega_{*e}$ ではこの近似で負であり、実部の微分 $\omega_{*e}/\omega_r^2>0$ に対して成長する。この局所模型の密度勾配駆動は :term[普遍的不安定性]{en="universal instability" yomi="ふへんてきふあんていせい"} の例である。磁気シア、有限イオン平行運動、温度勾配を含むと成長条件は変わる。

:::example[分極の効果]
$k_\perp\rho_s=1$ では $\omega_r=\omega_{*e}/2$ である。電子は $\omega-\omega_{*e}$ の向きに応答するので、この場合は電子の位相差が駆動になる。
:::

:::example[平行運動無視の条件]
イオン圧を無視しても、有限の $E_\parallel$ はイオンを加速できる。平行速度と連続式からその密度寄与は $\delta n_{i,\parallel}/n_0=(e\phi/T_e)(k_\parallel c_s/\omega)^2$ である。$k_\parallel c_s/\omega=0.1$ なら係数は $0.01$ なので、$\omega_{*e}/\omega$ が一程度のときにこれを小さいとできる。
:::

### 衝突による電子の位相差

同じスラブで、電子の慣性より衝突が強い極限を調べる。$\nu_{ei}\gg|\omega|$、$|k_\parallel|v_{te}\ll\nu_{ei}$ とし、平行運動でイオン流速を省く。$\Phi=e\phi/T_e$、$N=\delta n_e/n_0$ とすると、平行電気力と圧力勾配を摩擦に等置して $u_{e\parallel}=ik_\parallel T_e(\Phi-N)/(m_e\nu_{ei})$ である。連続式は $-i\omega N+i\omega_{*e}\Phi+ik_\parallel u_{e\parallel}=0$ なので、$D_\parallel=k_\parallel^2T_e/(m_e\nu_{ei})$ とおけば
$$
\frac{N}{\Phi}=\frac{D_\parallel-i\omega_{*e}}{D_\parallel-i\omega}
=1+\frac{i(\omega-\omega_{*e})}{D_\parallel}
+\order{D_\parallel^{-2}}
$$
となる。展開は $|\omega|,|\omega_{*e}|\ll D_\parallel$ の場合に用いる。冷イオンの密度と等置して、$k_\parallel c_s\ll|\omega|$ も保てば、実根は $\omega_r=\omega_{*e}/(1+k_\perp^2\rho_s^2)$、成長率は
$$
\gamma\simeq \frac{(\omega_{*e}-\omega_r)\omega_r^2}{D_\parallel\omega_{*e}} .
$$
これは衝突による :term[抵抗性ドリフト波]{en="resistive drift wave" yomi="ていこうせいどりふとは"} の位相差である。例えば $k_\perp\rho_s=1$ では $\gamma\simeq\omega_{*e}^2/(8D_\parallel)$ である。$\nu_{ei}$ をゼロへ近づけるとこの結果を使う条件が破れ、前の無衝突の速度積分が必要になる。 :cite[chen-plasma]

### 有限イオン分極でドリフト波が短波長化する

一様磁場、$x$ 方向の密度勾配、等温 Boltzmann 電子、冷イオンを仮定する。$E\times B$ の密度移流は $\omega_* e\phi/T_e$ に比例し、分極流の発散は $\omega k_\perp^2\rho_s^2e\phi/T_e$ に比例する。準中性により電子密度の時間微分 $\omega e\phi/T_e$ とつり合わせると $\omega(1+k_\perp^2\rho_s^2)=\omega_*$ となる。$k_x=0$ で $\omega_*=k_yV_*$ とすると群速は $V_*(1-k_y^2\rho_s^2)/(1+k_y^2\rho_s^2)^2$ になる。群速がゼロ、負になる場合でもこれは実周波数の変化であり、成長の証拠ではない。成長のためには粒子の応答の位相差が要る。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $x$ | $k_y\rho_s$ | 1 |
| $V_*$ | 反磁性角周波数を $k_y$ で割った値 | $\mathrm{m\,s^{-1}}$ |

分極込みの密度応答因子を $P$ とする。イオンの分極流を加える。

$$
P=x^{2} + 1
$$

規格化周波数を $\omega\rho_s/V_*$ とする。一次の反磁性周波数を分極で減らす。

$$
\omega\rho_s/V_*=\frac{x}{P}
$$

位相速度を $v_\phi$ とする。$x>0$ の例を考える。

$$
v_\phi=\frac{V_{*}}{P}
$$

群速度を $v_g$ とする。$k_y$ を変えて微分する。

$$
v_g=\frac{V_{*} \left(1 - x^{2}\right)}{P^{2}}
$$

二速度の符号を含む比を $v_g/v_\phi$ とする。短波長で負になる。

$$
v_g/v_\phi=\frac{v_{g}}{v_{\phi}}
$$

分極の割合を $f_{\rm pol}$ とする。応答因子のうち分極が担う割合である。

$$
f_{\rm pol}=\frac{x^{2}}{P}
$$

:::example[小さい波数]
反磁性速度と二つの速度を比較する。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{10} \\V_*&=100\,\mathrm{m\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

分極込みの密度応答因子は

$$
P=\left(\frac{1}{10}\right)^{2} + 1=\frac{101}{100}
$$

である。

規格化周波数は

$$
\omega\rho_s/V_*=\frac{\left(\frac{1}{10}\right)}{\left(\frac{101}{100}\right)}=\frac{10}{101}
$$

である。

位相速度は

$$
v_\phi=\frac{\left(100\right)}{\left(\frac{101}{100}\right)}=\frac{10000}{101}\,\mathrm{m\,s^{-1}}
$$

である。

群速度は

$$
v_g=\frac{\left(100\right) \left(1 - \left(\frac{1}{10}\right)^{2}\right)}{\left(\frac{101}{100}\right)^{2}}=\frac{990000}{10201}\,\mathrm{m\,s^{-1}}
$$

である。

二速度の符号を含む比は

$$
v_g/v_\phi=\frac{\left(\frac{990000}{10201}\right)}{\left(\frac{10000}{101}\right)}=\frac{99}{101}
$$

である。

分極の割合は

$$
f_{\rm pol}=\frac{\left(\frac{1}{10}\right)^{2}}{\left(\frac{101}{100}\right)}=\frac{1}{101}
$$

である。

分極補正は小さく、位相と群は同方向に進む。
:::

:::example[周波数の極大点]
群速度が消える波数を確かめる。

入力値を次のように置く。

$$
\begin{aligned}
x&=1 \\V_*&=100\,\mathrm{m\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

分極込みの密度応答因子は

$$
P=\left(1\right)^{2} + 1=2
$$

である。

規格化周波数は

$$
\omega\rho_s/V_*=\frac{\left(1\right)}{\left(2\right)}=\frac{1}{2}
$$

である。

位相速度は

$$
v_\phi=\frac{\left(100\right)}{\left(2\right)}=50\,\mathrm{m\,s^{-1}}
$$

である。

群速度は

$$
v_g=\frac{\left(100\right) \left(1 - \left(1\right)^{2}\right)}{\left(2\right)^{2}}=0\,\mathrm{m\,s^{-1}}
$$

である。

二速度の符号を含む比は

$$
v_g/v_\phi=\frac{\left(0\right)}{\left(50\right)}=0
$$

である。

分極の割合は

$$
f_{\rm pol}=\frac{\left(1\right)^{2}}{\left(2\right)}=\frac{1}{2}
$$

である。

位相は進んでも狭い波束の一次の群速度はゼロである。
:::

:::example[短波長側]
位相と包絡の向きが逆になるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=2 \\V_*&=100\,\mathrm{m\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

分極込みの密度応答因子は

$$
P=\left(2\right)^{2} + 1=5
$$

である。

規格化周波数は

$$
\omega\rho_s/V_*=\frac{\left(2\right)}{\left(5\right)}=\frac{2}{5}
$$

である。

位相速度は

$$
v_\phi=\frac{\left(100\right)}{\left(5\right)}=20\,\mathrm{m\,s^{-1}}
$$

である。

群速度は

$$
v_g=\frac{\left(100\right) \left(1 - \left(2\right)^{2}\right)}{\left(5\right)^{2}}=-12\,\mathrm{m\,s^{-1}}
$$

である。

二速度の符号を含む比は

$$
v_g/v_\phi=\frac{\left(-12\right)}{\left(20\right)}=- \frac{3}{5}
$$

である。

分極の割合は

$$
f_{\rm pol}=\frac{\left(2\right)^{2}}{\left(5\right)}=\frac{4}{5}
$$

である。

群速度の符号だけからエネルギーの増大は結論できない。
:::

:::example[反磁性の向きの反転]
密度勾配の反転の影響を調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=2 \\V_*&=-100\,\mathrm{m\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

分極込みの密度応答因子は

$$
P=\left(2\right)^{2} + 1=5
$$

である。

規格化周波数は

$$
\omega\rho_s/V_*=\frac{\left(2\right)}{\left(5\right)}=\frac{2}{5}
$$

である。

位相速度は

$$
v_\phi=\frac{\left(-100\right)}{\left(5\right)}=-20\,\mathrm{m\,s^{-1}}
$$

である。

群速度は

$$
v_g=\frac{\left(-100\right) \left(1 - \left(2\right)^{2}\right)}{\left(5\right)^{2}}=12\,\mathrm{m\,s^{-1}}
$$

である。

二速度の符号を含む比は

$$
v_g/v_\phi=\frac{\left(12\right)}{\left(-20\right)}=- \frac{3}{5}
$$

である。

分極の割合は

$$
f_{\rm pol}=\frac{\left(2\right)^{2}}{\left(5\right)}=\frac{4}{5}
$$

である。

位相と群の両方が反転し、相対的な向きは変化しない。
:::


## Weibel 不安定性

無磁場でも、粒子が一方向に多くの熱運動を持つと、横電流が磁場で集められる。これが :term[Weibel 不安定性]{en="Weibel instability" yomi="わいべるふあんていせい"} である。電子の $x,y$ 温度を $T_\perp$、$z$ 温度を $T_\parallel$ とし、固定イオンを置く。$k\hat{\vb z}$、$\vb E=E_x\hat{\vb x}$ なら $B_y=kE_x/\omega$。力と分布の速度勾配との積は
$$
(\vb E+\vb v\times\vb B_1)\cdot\nabla_v f_0
=-\frac{mE_xv_x}{T_\perp}f_0\left(1+A\frac{kv_z}{\omega}\right),
\quad A=T_\perp/T_\parallel-1.
$$
線形 Vlasov 式を $-i(\omega-kv_z)$ で割り、$qv_x$ のモーメントを取る。$\int v_x^2 f_0\,\dd v_x\dd v_y=(T_\perp/m)\int f_0\,\dd v_x\dd v_y$ と $\int F(v_z)/(\omega-kv_z)\dd v_z=-Z(\zeta)/(\sqrt2 kv_{t\parallel})$ を使う。Ampère–Maxwell 則へ入れると
$$
\omega^2-c^2k^2+\omega_p^2[A+(1+A)\zeta Z(\zeta)]=0.
\label{eq-plasma-weibel}
$$
$\omega=i\gamma,\ \gamma>0$ なら $\zeta=iy$ であり、$\zeta Z(\zeta)=-\int e^{-x^2}y^2/(x^2+y^2)\dd x/\sqrt\pi$ である。これは $y=0$ でゼロ、増加する $y$ で単調減少して $-1$ に近づく。左辺は $\gamma$ とともに単調減少するため、$0<c^2k^2<\omega_p^2 A$ なら純成長根が一個ある。逆にこの温度分布には、同じ伝播方向に純成長以外の不安定根はない。実際、$\Im\zeta>0$ で偶関数の重みを用いると $\zeta Z(\zeta)$ は $\zeta^2/(x^2-\zeta^2)$ の Maxwell 平均になる。この被積分関数の虚部は $x^2\Im(\zeta^2)/|x^2-\zeta^2|^2$ である。したがって $\Im(\omega^2)\ne0$ なら $\omega^2$ とこの平均の虚部は同符号であり、分散関係を満たせない。$\Im\omega>0$ のもとで $\Im(\omega^2)=0$ なら $\Re\omega=0$ である。

:::example[不安定な波数]
$A=1/4$ では $0<kd_e<1/2$ が純成長分枝の範囲である。$kd_e=1/2$ は境界である。
:::

### しきい値の近傍

$\zeta=iy$、$y=\gamma/(\sqrt2 k v_{t\parallel})\ll1$ なら $\zeta Z(\zeta)=-\sqrt\pi y+\order{y^2}$ である。式 :eqref[eq-plasma-weibel] で $\gamma^2$ も高次として落とすと
$$
\gamma\simeq\sqrt{\frac2\pi}k v_{t\parallel}\frac{A-k^2d_e^2}{1+A}.
$$
これはしきい値近くの正の差で使える。温度異方性がなくなるまでの非線形発展は、この平衡分布を固定した式からは求まらない。

:::example[限界に近い波数]
$A=1/4$、$k^2d_e^2=0.24$ なら $\gamma/(k v_{t\parallel})\simeq0.008\sqrt{2/\pi}$ である。この値の小ささは $y\ll1$ と整合する。
:::

### Weibel の限界近傍の漸近幅

一方向に伝播する横波で、波の電場方向の温度と波数方向の温度の比を $1+A$ とする。分散式は $D=\omega^2-c^2k^2+\omega_p^2[A+(1+A)\zeta Z(\zeta)]$ である。$\omega=i\gamma$ とし $\zeta=iy$、$y=\gamma/(\sqrt2 k v_t)$ とおくと $\zeta Z=-\sqrt\pi y e^{y^2}\operatorname{erfc}y$ である。$y\ll1$ かつ $\gamma^2/\omega_p^2\ll A-(kd_e)^2$ なら $\gamma\simeq\sqrt{2/\pi}\,k v_t[A-(kd_e)^2]/(1+A)$ を得る。$0<(kd_e)^2<A$ で正になる。近似の成長率の微分から $kd_e=\sqrt{A/3}$ で最大だが、$A$ が大きい場合にはこの位置を全式の最大点とみなさない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $A$ | 温度異方性の過剰分 | 1 |
| $x$ | 波数 $kd_e$ | 1 |
| $U$ | 伝播方向の熱速さ $v_t/c$ | 1 |

しきい値からの距離を $H$ とする。正の側で非振動的な成長を調べる。

$$
H=A - x^{2}
$$

小引数の虚周波数を $y$ とする。$Z(0)=i\sqrt\pi$ を使う。

$$
y=\frac{H}{\sqrt{\pi} \left(A + 1\right)}
$$

規格化した成長率の近似を $\gamma/\omega_p$ とする。負なら不安定枝の推定として使わない。

$$
\gamma/\omega_p=\sqrt{2} U x y
$$

近似の最大点を $k_* d_e$ とする。波数の三次式を微分する。

$$
k_* d_e=\frac{\sqrt{3} \sqrt{A}}{3}
$$

近似の最大成長率を $G_*$ とする。最大波数で評価する。

$$
G_*=\frac{2\sqrt2 A U(k_*d_e)}{3\sqrt\pi(1+A)}
$$

捨てた慣性項の相対尺度を $\mathcal E$ とする。$H\ne0$ で比較する。

$$
\mathcal E=\frac{(\gamma/\omega_p)^2}{|H|}
$$

:::example[しきい値より十分下の低波数]
駆動の大きさと成長率を比較する。

入力値を次のように置く。

$$
\begin{aligned}
A&=\frac{1}{10} \\x&=\frac{1}{20} \\U&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

しきい値からの距離は

$$
H=\left(\frac{1}{10}\right) - \left(\frac{1}{20}\right)^{2}=\frac{39}{400}
$$

である。

小引数の虚周波数は

$$
y=\frac{\left(\frac{39}{400}\right)}{\sqrt{\pi} \left(\left(\frac{1}{10}\right) + 1\right)}\simeq 0.050008
$$

である。

規格化した成長率の近似は

$$
\gamma/\omega_p=\sqrt{2} \left(\frac{1}{100}\right) \left(\frac{1}{20}\right) \left(\frac{39}{440 \sqrt{\pi}}\right)\simeq 3.5361\times 10^{-5}
$$

である。

近似の最大点は

$$
k_* d_e=\frac{\sqrt{3} \sqrt{\left(\frac{1}{10}\right)}}{3}\simeq 0.18257
$$

である。

近似の最大成長率は

$$
G_*=\frac{2 \sqrt{2} \left(\frac{1}{100}\right) \left(\frac{1}{10}\right) \left(\frac{\sqrt{30}}{30}\right)}{3 \sqrt{\pi} \left(\left(\frac{1}{10}\right) + 1\right)}\simeq 8.8287\times 10^{-5}
$$

である。

捨てた慣性項の相対尺度は

$$
\mathcal E=\frac{\left(\frac{39 \sqrt{2}}{880000 \sqrt{\pi}}\right)^{2}}{\left|{\left(\frac{39}{400}\right)}\right|}\simeq 1.2824\times 10^{-8}
$$

である。

駆動が有限でも $k v_t$ が小さいため成長は遅い。
:::

:::example[近似の最速位置]
最大値との一致を求める。

入力値を次のように置く。

$$
\begin{aligned}
A&=\frac{3}{100} \\x&=\frac{1}{10} \\U&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

しきい値からの距離は

$$
H=- \left(\frac{1}{10}\right)^{2} + \left(\frac{3}{100}\right)=\frac{1}{50}
$$

である。

小引数の虚周波数は

$$
y=\frac{\left(\frac{1}{50}\right)}{\sqrt{\pi} \left(\left(\frac{3}{100}\right) + 1\right)}\simeq 0.010955
$$

である。

規格化した成長率の近似は

$$
\gamma/\omega_p=\sqrt{2} \left(\frac{1}{100}\right) \left(\frac{1}{10}\right) \left(\frac{2}{103 \sqrt{\pi}}\right)\simeq 1.5493\times 10^{-5}
$$

である。

近似の最大点は

$$
k_* d_e=\frac{\sqrt{3} \sqrt{\left(\frac{3}{100}\right)}}{3}=\frac{1}{10}
$$

である。

近似の最大成長率は

$$
G_*=\frac{2 \sqrt{2} \left(\frac{1}{100}\right) \left(\frac{1}{10}\right) \left(\frac{3}{100}\right)}{3 \sqrt{\pi} \left(\left(\frac{3}{100}\right) + 1\right)}\simeq 1.5493\times 10^{-5}
$$

である。

捨てた慣性項の相対尺度は

$$
\mathcal E=\frac{\left(\frac{\sqrt{2}}{51500 \sqrt{\pi}}\right)^{2}}{\left|{\left(\frac{1}{50}\right)}\right|}\simeq 1.2002\times 10^{-8}
$$

である。

選んだ波数は近似式の最速波数である。
:::

:::example[短波長の安定側]
形式の成長率と入力の条件を照合する。

入力値を次のように置く。

$$
\begin{aligned}
A&=\frac{3}{100} \\x&=\frac{1}{5} \\U&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

しきい値からの距離は

$$
H=- \left(\frac{1}{5}\right)^{2} + \left(\frac{3}{100}\right)=- \frac{1}{100}
$$

である。

小引数の虚周波数は

$$
y=\frac{\left(- \frac{1}{100}\right)}{\sqrt{\pi} \left(\left(\frac{3}{100}\right) + 1\right)}\simeq -0.0054776
$$

である。

規格化した成長率の近似は

$$
\gamma/\omega_p=\sqrt{2} \left(- \frac{1}{103 \sqrt{\pi}}\right) \left(\frac{1}{100}\right) \left(\frac{1}{5}\right)\simeq -1.5493\times 10^{-5}
$$

である。

近似の最大点は

$$
k_* d_e=\frac{\sqrt{3} \sqrt{\left(\frac{3}{100}\right)}}{3}=\frac{1}{10}
$$

である。

近似の最大成長率は

$$
G_*=\frac{2 \sqrt{2} \left(\frac{1}{100}\right) \left(\frac{1}{10}\right) \left(\frac{3}{100}\right)}{3 \sqrt{\pi} \left(\left(\frac{3}{100}\right) + 1\right)}\simeq 1.5493\times 10^{-5}
$$

である。

捨てた慣性項の相対尺度は

$$
\mathcal E=\frac{\left(- \frac{\sqrt{2}}{51500 \sqrt{\pi}}\right)^{2}}{\left|{\left(- \frac{1}{100}\right)}\right|}\simeq 2.4003\times 10^{-8}
$$

である。

駆動が負なので正成長の近似として解釈できない。
:::

:::example[温度方向の偏りが大きい場合]
小引数展開に十分な分離があるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
A&=10 \\x&=1 \\U&=\frac{1}{10}
\end{aligned}
$$

これらを順に代入する。

しきい値からの距離は

$$
H=\left(10\right) - \left(1\right)^{2}=9
$$

である。

小引数の虚周波数は

$$
y=\frac{\left(9\right)}{\sqrt{\pi} \left(\left(10\right) + 1\right)}\simeq 0.46161
$$

である。

規格化した成長率の近似は

$$
\gamma/\omega_p=\sqrt{2} \left(1\right) \left(\frac{1}{10}\right) \left(\frac{9}{11 \sqrt{\pi}}\right)\simeq 0.065281
$$

である。

近似の最大点は

$$
k_* d_e=\frac{\sqrt{3} \sqrt{\left(10\right)}}{3}\simeq 1.8257
$$

である。

近似の最大成長率は

$$
G_*=\frac{2 \sqrt{2} \left(10\right) \left(\frac{1}{10}\right) \left(\frac{\sqrt{30}}{3}\right)}{3 \sqrt{\pi} \left(\left(10\right) + 1\right)}\simeq 0.088287
$$

である。

捨てた慣性項の相対尺度は

$$
\mathcal E=\frac{\left(\frac{9 \sqrt{2}}{110 \sqrt{\pi}}\right)^{2}}{\left|{\left(9\right)}\right|}\simeq 0.00047352
$$

である。

虚引数が小さいとは言いにくく、誤差尺度に低次の慣性項だけを使うと不十分である。
:::


## 温度異方性による不安定性

背景磁場がある場合、圧力異方性は磁気張力を変える。二重断熱の圧力を $\mathsf P=p_\perp I+(p_\parallel-p_\perp)\vb b\vb b$ とする。一様平衡の平行横波で $\delta\vb b_\perp=-k_\parallel B_0\vb u_\perp/(\omega B_0)$ である。圧力力の横成分は $-ik_\parallel(p_\parallel-p_\perp)\delta\vb b_\perp$ であり、磁気張力と足すと
$$
\omega^2=\frac{k_\parallel^2}{\rho}
\left(\frac{B_0^2}{\mu_0}+p_\perp-p_\parallel\right).
\label{eq-plasma-firehose}
$$
$p_\parallel-p_\perp>B_0^2/\mu_0$ で :term[ファイアホース不安定性]{en="firehose instability" yomi="ふぁいあほーすふあんていせい"} が起こる。この式は長波長の平行流体模型のしきい値である。

$p_\perp>p_\parallel$ では :term[ミラー不安定性]{en="mirror instability" yomi="みらーふあんていせい"} という圧縮性の不安定性もある。一つの異方イオン種と冷たい中和電子を仮定し、低周波、長波長、$0<|k_\parallel|\ll k_\perp$ の運動論極限ではしきい値は $\beta_\perp(T_\perp/T_\parallel-1)>1$ である。ただし $\beta_\perp=2\mu_0p_\perp/B_0^2$ とする。鏡力で粒子が弱い磁場へ集まり、垂直圧力が磁気圧の変化を補強する。一般の電子応答、有限旋回半径、斜めファイアホースのしきい値は異方性の運動論固有値問題が必要であり、ここでは平行流体分枝を導出し、鏡型の扱いは :cite[stix-waves] を参照する。

:::example[張力の反転]
$p_\parallel-p_\perp=2B_0^2/\mu_0$ なら式 :eqref[eq-plasma-firehose] は $\omega^2=-k_\parallel^2 v_A^2$ である。正の時間成長率は $\abs{k_\parallel}v_A$ になる。
:::

:::example[平行ベータでの表現]
$\beta_\parallel=2\mu_0p_\parallel/B_0^2$ を定義すると、長波長ファイアホース条件は $1-p_\perp/p_\parallel>2/\beta_\parallel$ である。$\beta_\parallel=8$ なら $p_\perp/p_\parallel<3/4$ が条件となる。高い熱圧では比較的小さい比の異方性でも張力が打ち消されうる。
:::

### 鏡型のしきい値の静的導出

低周波のドリフト運動論で、局所一様な磁場を $\vb B_0=B_0\hat{\vb z}$、冷たい質量なしの電子と二温度 Maxwell イオンを仮定する。電子圧を落とすため $E_\parallel=0$ である。$0<\abs{k_\parallel}\ll k_\perp$ とし、臨界で $\omega/(k_\parallel v_{ti\parallel})\to0$ を先にとる。磁気モーメントを固定したイオン分布は $F_0\propto\exp[-m_i v_\parallel^2/(2T_\parallel)-\mu B_0/T_\perp]$ である。線形ドリフト運動論の平行移流と鏡力は
$$
v_\parallel i k_\parallel\,\delta F
-\frac{\mu}{m_i}i k_\parallel\delta B\,\partial_{v_\parallel}F_0=0
$$
であり、連続な静的応答は $\delta F=-\mu\delta B\,F_0/T_\parallel$ である。速度積分の測度 $2\pi B\,\dd\mu\,\dd v_\parallel/m_i$ も摂動する。垂直圧力は $\int\mu B\, F\,\dd^3v$ なので、$B$ の二因子から $2p_\perp\delta B/B_0$ が出る。分布の補正は $\mu$ の指数積分で $\langle\mu^2\rangle/\langle\mu\rangle=2T_\perp/B_0$ を使うと $-2p_\perp(T_\perp/T_\parallel)\delta B/B_0$ である。従って
$$
\delta\left(p_\perp+\frac{B^2}{2\mu_0}\right)
=\left[\frac{B_0^2}{\mu_0}-2p_\perp\left(\frac{T_\perp}{T_\parallel}-1\right)\right]\frac{\delta B}{B_0}.
$$
$k_\parallel/k_\perp\to0$ で、垂直力の主要項はこの全圧の勾配である。括弧が負に変わると、弱い場が粒子圧の増加によってさらに弱くなる応答を持つ。この極限の中立条件から $\beta_\perp(T_\perp/T_\parallel-1)=1$ を得る。有限の傾きでの磁気張力と有限旋回半径は、成長率と波数範囲を変える。そこでこれらを落としたまま平行の速度積分だけを残し、中立条件のどちら側に成長根があるかを確かめる。$k_\parallel>0$、$\zeta=\omega/(\sqrt2 k_\parallel v_{ti\parallel})$、$A_i=T_\perp/T_\parallel$ とする。平行移流を $v_\parallel-\omega/k_\parallel$ に置き換えると、$\delta F$ の静的な式に $v_\parallel/(v_\parallel-\omega/k_\parallel)$ が掛かる。その Maxwell 平均は $1+\zeta Z(\zeta)$ なので、垂直全圧の釣り合いは
$$
1+\beta_\perp\{1-A_i[1+\zeta Z(\zeta)]\}=0
$$
となる。$\zeta=iy$、$y>0$ では $1+\zeta Z(\zeta)$ が $1$ から $0$ へ単調に下がることを Weibel の節の積分表示から確認できる。従って静的なしきい値を超えた場合にこの式は純成長根を一つ持つ。しきい値の近くで $0<\beta_\perp(A_i-1)-1\ll \beta_\perp A_i$ のとき
$$
\gamma\simeq
\sqrt{\frac{2}{\pi}}\,k_\parallel v_{ti\parallel}
\frac{\beta_\perp(A_i-1)-1}{\beta_\perp A_i}
$$
である。この式が与えるのは、慣性・張力・有限旋回半径を落とした低周波の先頭項である。

:::example[磁化の強さを比べる]
同じ $T_\perp/T_\parallel=1.2$ でも、$\beta_\perp=4$ なら積は $0.8$、$\beta_\perp=6$ なら $1.2$ である。後者はこの低周波模型で鏡型のしきい値を越える。$p_\parallel<p_\perp$ なので、平行ファイアホースの条件とは同時に成り立たない。
:::

### 異方性張力で平行波が反転する条件

ジャイロ対称圧力 $\mathsf P=p_\perp I+(p_\parallel-p_\perp)\hat{\vb b}\hat{\vb b}$ に対する磁力線の復元係数は $B^2/\mu_0+p_\perp-p_\parallel$ である。実際、横変位 $\xi_\perp(z)$ の一次で $\delta\hat{\vb b}=\partial_z\xi_\perp$ となり、$-\divg\mathsf P$ は $-(p_\parallel-p_\perp)\partial_z^2\xi_\perp$ を与える。$\vb J\times\vb B$ の張力項と合わせて $\omega^2/k_\parallel^2=v_A^2[1+\beta_\parallel(A-1)/2]$ を得る。ここで $A=p_\perp/p_\parallel$ である。平行 Firehose 条件だけを与える長波長流体式であり、斜めの運動論的不安定性や実際の限界波数は別に求める。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $A$ | 平行圧に対する垂直圧の比 | 1 |
| $\beta_\parallel$ | 平行ベータ | 1 |
| $K$ | 平行波数と基準長の積 | 1 |

垂直ベータを $\beta_\perp$ とする。磁場強度の単位は共通とする。

$$
\beta_\perp=A \beta_{\parallel}
$$

規格化張力を $\mathcal T$ とする。圧力異方性の方向で増減する。

$$
\mathcal T=\frac{\beta_{\parallel} \left(A - 1\right)}{2} + 1
$$

張力がゼロになる異方性を $A_c$ とする。$\beta_\parallel>0$ で求める。

$$
A_c=1 - \frac{2}{\beta_{\parallel}}
$$

しきい値との差を $A-A_c$ とする。負なら張力が負である。

$$
A-A_c=A - A_{c}
$$

規格化した周波数二乗を $W^2$ とする。単位は Alfvén 速度を基準長で割った値の二乗である。

$$
W^2=K^{2} \mathcal T
$$

単一異方種の単純な鏡型指標を $\mathcal M$ とする。長波長の鏡型条件との比較用で、Firehose とは独立である。

$$
\mathcal M=\beta_{\perp} \left(A - 1\right) - 1
$$

:::example[等方平衡]
Alfvén 波の張力が戻るか確かめる。

入力値を次のように置く。

$$
\begin{aligned}
A&=1 \\\beta_\parallel&=10 \\K&=1
\end{aligned}
$$

これらを順に代入する。

垂直ベータは

$$
\beta_\perp=\left(10\right) \left(1\right)=10
$$

である。

規格化張力は

$$
\mathcal T=\frac{\left(10\right) \left(\left(1\right) - 1\right)}{2} + 1=1
$$

である。

張力がゼロになる異方性は

$$
A_c=1 - \frac{2}{\left(10\right)}=\frac{4}{5}
$$

である。

しきい値との差は

$$
A-A_c=\left(1\right) - \left(\frac{4}{5}\right)=\frac{1}{5}
$$

である。

規格化した周波数二乗は

$$
W^2=\left(1\right)^{3}=1
$$

である。

単一異方種の単純な鏡型指標は

$$
\mathcal M=\left(10\right) \left(\left(1\right) - 1\right) - 1=-1
$$

である。

異方性補正は消える。
:::

:::example[張力が消える場合]
高ベータの小さな異方性の効果を求める。

入力値を次のように置く。

$$
\begin{aligned}
A&=\frac{4}{5} \\\beta_\parallel&=10 \\K&=1
\end{aligned}
$$

これらを順に代入する。

垂直ベータは

$$
\beta_\perp=\left(10\right) \left(\frac{4}{5}\right)=8
$$

である。

規格化張力は

$$
\mathcal T=\frac{\left(10\right) \left(\left(\frac{4}{5}\right) - 1\right)}{2} + 1=0
$$

である。

張力がゼロになる異方性は

$$
A_c=1 - \frac{2}{\left(10\right)}=\frac{4}{5}
$$

である。

しきい値との差は

$$
A-A_c=0
$$

である。

規格化した周波数二乗は

$$
W^2=\left(0\right) \left(1\right)^{2}=0
$$

である。

単一異方種の単純な鏡型指標は

$$
\mathcal M=\left(8\right) \left(\left(\frac{4}{5}\right) - 1\right) - 1=- \frac{13}{5}
$$

である。

磁気張力と異方圧力が相殺する。
:::

:::example[平行圧が優勢]
負の二乗周波数を確認する。

入力値を次のように置く。

$$
\begin{aligned}
A&=\frac{1}{2} \\\beta_\parallel&=10 \\K&=1
\end{aligned}
$$

これらを順に代入する。

垂直ベータは

$$
\beta_\perp=\left(10\right) \left(\frac{1}{2}\right)=5
$$

である。

規格化張力は

$$
\mathcal T=\frac{\left(10\right) \left(\left(\frac{1}{2}\right) - 1\right)}{2} + 1=- \frac{3}{2}
$$

である。

張力がゼロになる異方性は

$$
A_c=1 - \frac{2}{\left(10\right)}=\frac{4}{5}
$$

である。

しきい値との差は

$$
A-A_c=\left(\frac{1}{2}\right) - \left(\frac{4}{5}\right)=- \frac{3}{10}
$$

である。

規格化した周波数二乗は

$$
W^2=\left(- \frac{3}{2}\right) \left(1\right)^{2}=- \frac{3}{2}
$$

である。

単一異方種の単純な鏡型指標は

$$
\mathcal M=\left(5\right) \left(\left(\frac{1}{2}\right) - 1\right) - 1=- \frac{7}{2}
$$

である。

この長波長の平行分枝は指数成長側である。
:::

:::example[垂直圧が優勢]
平行張力と鏡型指標を比較する。

入力値を次のように置く。

$$
\begin{aligned}
A&=2 \\\beta_\parallel&=10 \\K&=1
\end{aligned}
$$

これらを順に代入する。

垂直ベータは

$$
\beta_\perp=\left(10\right) \left(2\right)=20
$$

である。

規格化張力は

$$
\mathcal T=\frac{\left(10\right) \left(\left(2\right) - 1\right)}{2} + 1=6
$$

である。

張力がゼロになる異方性は

$$
A_c=1 - \frac{2}{\left(10\right)}=\frac{4}{5}
$$

である。

しきい値との差は

$$
A-A_c=\left(2\right) - \left(\frac{4}{5}\right)=\frac{6}{5}
$$

である。

規格化した周波数二乗は

$$
W^2=\left(1\right)^{2} \left(6\right)=6
$$

である。

単一異方種の単純な鏡型指標は

$$
\mathcal M=\left(20\right) \left(\left(2\right) - 1\right) - 1=19
$$

である。

平行 Firehose の復元は強くなるが、異なる偏波の鏡型条件は満たし得る。
:::


## 演習問題

::::exercise{level=basic}
冷たい対称ビームで $kU=2\omega_p$ のとき負の $X$ はあるか。

:::solution
$a=4\omega_p^2$ で積は正、和は正、判別式も正なので二つとも正である。負の $X$ はない。
:::
::::

::::exercise{level=basic}
一峰の Maxwell 分布の極大で $P_F$ の符号を求めよ。

:::solution
全速度で $F(v)\le F(u_0)$ なので積分は負である。
:::
::::

::::exercise{level=basic}
$k_\perp\rho_s=2$ のドリフト波周波数を求めよ。

:::solution
$\omega/\omega_{*e}=1/(1+4)=1/5$ である。
:::
::::

::::exercise{level=basic}
$A=1$ での Weibel の上限 $kd_e$ を求めよ。

:::solution
しきい値 $k^2d_e^2=A$ より $kd_e=1$ である。
:::
::::

::::exercise{level=advanced}
対称ビームの最大成長率を微分で確かめよ。

:::solution
$\dv{X_-}{a}=1-2\omega_p^2/\sqrt{\omega_p^4+8a\omega_p^2}$ で、ゼロは $a=3\omega_p^2/8$ である。微分が負から正に変わる。これを $X_-$ に戻すと $-\omega_p^2/8$ である。
:::
::::

::::exercise{level=advanced}
冷たいビーム密度を八倍にしたとき希薄ビーム成長率を比較せよ。

:::solution
$\omega_b^2\propto n_b$、共鳴条件の成長率は $(\omega_b^2)^{1/3}$ に比例するため二倍である。背景の密度を固定し、希薄条件を保つ。
:::
::::

::::exercise{level=advanced}
異方性がゼロなら Weibel 式の純成長根は存在するか。

:::solution
$A=0$ では $\gamma=0$ で左辺は $-c^2k^2<0$ である。$\gamma>0$ で各追加項は非正なのでゼロには到達しない。
:::
::::

::::exercise{level=advanced}
$p_\parallel-p_\perp=B_0^2/(2\mu_0)$ で平行波の位相速度を求めよ。

:::solution
残った張力は元の半分である。式 :eqref[eq-plasma-firehose] の正の位相速度は $v_A/\sqrt2$ となる。
:::
::::

## 確認テスト

::::quiz
同じ電子温度なら、密度勾配を強くして局所ドリフト波はどう変わるか。

- [x] $L_n^{-1}$ に比例して角周波数が増す。
- [ ] 必ず周波数がゼロになる。
- [ ] 周波数は勾配に依存しない。

:::explanation
固定波数で $\omega_{*e}=k_y T_e/(eBL_n)$、分極因子は $L_n$ に依存しない。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$A=4$ の Weibel の限界波数 $kd_e$ を答えよ。

:::explanation
$\omega_p^2 A=c^2 k^2$ である。正の $kd_e$ は $\sqrt4=2$ となる。
:::
::::

::::quiz
長波長平行ファイアホースを駆動する圧力差を選べ。

- [x] 十分大きい $p_\parallel-p_\perp>0$
- [ ] 十分大きい $p_\perp-p_\parallel>0$
- [ ] 等方圧力

:::explanation
圧力差が磁気張力 $B_0^2/\mu_0$ を上回ると式 :eqref[eq-plasma-firehose] の $\omega^2$ が負になる。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
$k_\perp\rho_s=1$ で $\omega/\omega_{*e}$ を答えよ。

:::explanation
断熱電子で比は $1/(1+k_\perp^2\rho_s^2)=1/2$ である。
:::
::::
