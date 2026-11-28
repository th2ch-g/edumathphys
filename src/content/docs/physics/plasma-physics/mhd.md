---
title: 電磁流体力学
description: 一流体の誘導と運動量から、MHD 波とエネルギー収支を導く。
notation:
  - symbol: '\rho'
    meaning: 質量密度
  - symbol: '\eta_m'
    meaning: 磁気拡散率
  - symbol: 'v_A'
    meaning: Alfvén 速度
  - symbol: '\vb\xi'
    meaning: 流体変位
---

## この章の目標

- 理想 MHD の仮定を明示する。
- 磁束保存と磁気応力を導く。
- 波と安定性をエネルギーから調べる。

## 理想 MHD 方程式

粒子尺度より十分に長い長さで、質量流と磁場を残す。前章までの二流体の和と :eqref[eq-plasma-ohm] で、準中性、変位電流無視、等方圧と断熱閉包を用いる。Hall 項、電子圧力による誘導、電子慣性、抵抗を無視すると $\vb E=-\vb U\times\vb B$ である。よって :term[理想電磁流体力学]{en="ideal magnetohydrodynamics" yomi="りそうでんじりゅうたいりきがく"} は
$$
\begin{aligned}
&\partial_t\rho+\nabla\cdot(\rho\vb U)=0,\quad
\rho D_t\vb U=-\nabla p+\vb j\times\vb B,\quad
D_t(p\rho^{-\gamma})=0,\\
&\partial_t\vb B=\nabla\times(\vb U\times\vb B),\quad
\nabla\cdot\vb B=0,\quad \mu_0\vb j=\nabla\times\vb B
\end{aligned}
\label{eq-plasma-ideal-mhd}
$$
で閉じる。$D_t=\partial_t+\vb U\cdot\nabla$ である。磁場の発散の時間微分はゼロであり、初期条件で発散をゼロにする。

### 総エネルギー

圧力式を内部エネルギー $p/(\gamma-1)$ の式に直す。運動式と $\vb U$ の内積では磁気力の仕事が $\vb U\cdot(\vb j\times\vb B)=\vb j\cdot\vb E$ である。これは :eqref[eq-em-poynting] の磁気部分の消失率と相殺する。従って $\mathcal E=\rho U^2/2+p/(\gamma-1)+B^2/(2\mu_0)$ は、流束 $(\rho U^2/2+\gamma p/(\gamma-1))\vb U+\vb E\times\vb B/\mu_0$ で局所保存される。

:::example[一様静止平衡]
$\rho,p,\vb B$ が一定で $\vb U=0$ なら $\vb j=0$、圧力勾配もゼロである。すべての時間微分が消え、波の基準平衡に使える。
:::

:::example[磁場に沿う定常せん断流]
$\rho$ 一定で $\vb U=U(z)\hat{\vb x}$、$\vb B=B(z)\hat{\vb x}$ とする。移流加速度と電場はゼロである。磁場に垂直な全圧を一定にする $p(z)=p_{\rm tot}-B(z)^2/(2\mu_0)$ を選べば定常な力の釣り合いが成立する。密度一定でも $p$ は流線を横切って変化でき、断熱式は $U\partial_xp=0$ を満たす。
:::

### 横圧縮を支える二種類の圧力

等方で断熱な一様平衡に対し、横にだけ圧縮する理想流体要素を取る。連続式と誘導式から $\dd\ln\rho/\dd t=\dd\ln B/\dd t$ である。断熱式は $p=p_0(\rho/\rho_0)^\gamma$ なので、体積要素の全圧力は $p_0 N^\gamma+B_0^2N^2/(2\mu_0)$ となる。$\rho$ で微分した全圧力の変化は $c_s^2+v_A^2$ であり、垂直の速い磁気音波の速度に一致する。磁場の強度や温度を一定にして微分するとこの復元力を得られない。横圧縮の条件では熱圧力と磁気圧力の応答が異なる指数で増える。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $N$ | 密度圧縮率 | 1 |
| $\beta_0$ | 初期の熱圧力と磁気圧力の比 | 1 |
| $\gamma$ | 断熱指数 | 1 |

熱圧力比を $p/p_0$ とする。$p/\rho^\gamma$ は一定である。

$$
p/p_0=N^{\gamma}
$$

磁気圧力比を $p_B/p_{B0}$ とする。磁束凍結により横圧縮で $B$ が増す。

$$
p_B/p_{B0}=N^{2}
$$

全圧力の値を $(p+p_B)/p_{B0}$ とする。熱圧力を磁気圧力単位へ換算する。

$$
(p+p_B)/p_{B0}=\beta_{0} p/p_{0} + p_B/p_{B0}
$$

圧縮後のベータを $\beta$ とする。式では圧力の比を用いる。

$$
\beta=\frac{\beta_{0} p/p_{0}}{p_B/p_{B0}}
$$

音速の二乗を $c_s^2/v_{A0}^2$ とする。$c_s^2=\gamma p/\rho$ である。

$$
c_s^2/v_{A0}^2=\frac{\beta_{0} \gamma p/p_{0}}{2 N}
$$

垂直速波の速度の二乗を $c_f^2/v_{A0}^2$ とする。$v_A^2/v_{A0}^2=N$ が加わる。

$$
c_f^2/v_{A0}^2=N + c_s^2/v_{A0}^2
$$

:::example[磁気圧が主な圧縮]
圧縮で熱圧が支配側へ移るか求める。

入力値を次のように置く。

$$
\begin{aligned}
N&=2 \\\beta_0&=\frac{1}{10} \\\gamma&=\frac{5}{3}
\end{aligned}
$$

これらを順に代入する。

熱圧力比は

$$
p/p_0=\left(2\right)^{\left(\frac{5}{3}\right)}\simeq 3.1748
$$

である。

磁気圧力比は

$$
p_B/p_{B0}=\left(2\right)^{2}=4
$$

である。

全圧力の値は

$$
(p+p_B)/p_{B0}=\left(2 \cdot 2^{\frac{2}{3}}\right) \left(\frac{1}{10}\right) + \left(4\right)\simeq 4.3175
$$

である。

圧縮後のベータは

$$
\beta=\frac{\left(2 \cdot 2^{\frac{2}{3}}\right) \left(\frac{1}{10}\right)}{\left(4\right)}\simeq 0.07937
$$

である。

音速の二乗は

$$
c_s^2/v_{A0}^2=\frac{\left(2 \cdot 2^{\frac{2}{3}}\right) \left(\frac{1}{10}\right) \left(\frac{5}{3}\right)}{2 \left(2\right)}\simeq 0.13228
$$

である。

垂直速波の速度の二乗は

$$
c_f^2/v_{A0}^2=\left(2\right) + \left(\frac{2^{\frac{2}{3}}}{12}\right)\simeq 2.1323
$$

である。

磁気圧の伸びの方が大きいのでベータはさらに小さくなる。
:::

:::example[熱圧が主な圧縮]
同じ密度変化で復元速度を比較する。

入力値を次のように置く。

$$
\begin{aligned}
N&=2 \\\beta_0&=10 \\\gamma&=\frac{5}{3}
\end{aligned}
$$

これらを順に代入する。

熱圧力比は

$$
p/p_0=\left(2\right)^{\left(\frac{5}{3}\right)}\simeq 3.1748
$$

である。

磁気圧力比は

$$
p_B/p_{B0}=\left(2\right)^{2}=4
$$

である。

全圧力の値は

$$
(p+p_B)/p_{B0}=\left(10\right) \left(2 \cdot 2^{\frac{2}{3}}\right) + \left(4\right)\simeq 35.748
$$

である。

圧縮後のベータは

$$
\beta=\frac{\left(10\right) \left(2 \cdot 2^{\frac{2}{3}}\right)}{\left(4\right)}\simeq 7.937
$$

である。

音速の二乗は

$$
c_s^2/v_{A0}^2=\frac{\left(10\right) \left(2 \cdot 2^{\frac{2}{3}}\right) \left(\frac{5}{3}\right)}{2 \left(2\right)}\simeq 13.228
$$

である。

垂直速波の速度の二乗は

$$
c_f^2/v_{A0}^2=\left(2\right) + \left(\frac{25 \cdot 2^{\frac{2}{3}}}{3}\right)\simeq 15.228
$$

である。

圧縮後も熱圧の復元が大きい。
:::

:::example[二つの指数が等しい模型]
圧縮してもベータが変わらない条件を確かめる。

入力値を次のように置く。

$$
\begin{aligned}
N&=2 \\\beta_0&=1 \\\gamma&=2
\end{aligned}
$$

これらを順に代入する。

熱圧力比は

$$
p/p_0=\left(2\right)^{\left(2\right)}=4
$$

である。

磁気圧力比は

$$
p_B/p_{B0}=\left(2\right)^{2}=4
$$

である。

全圧力の値は

$$
(p+p_B)/p_{B0}=\left(1\right) \left(4\right) + \left(4\right)=8
$$

である。

圧縮後のベータは

$$
\beta=\left(1\right)=1
$$

である。

音速の二乗は

$$
c_s^2/v_{A0}^2=\frac{\left(1\right) \left(4\right)}{2}=2
$$

である。

垂直速波の速度の二乗は

$$
c_f^2/v_{A0}^2=2 \left(2\right)=4
$$

である。

二つの圧力が同じ指数で増える。
:::

:::example[横への希釈]
磁場が薄まる場合に圧力比を調べる。

入力値を次のように置く。

$$
\begin{aligned}
N&=\frac{1}{8} \\\beta_0&=1 \\\gamma&=\frac{5}{3}
\end{aligned}
$$

これらを順に代入する。

熱圧力比は

$$
p/p_0=\left(\frac{1}{8}\right)^{\left(\frac{5}{3}\right)}=\frac{1}{32}
$$

である。

磁気圧力比は

$$
p_B/p_{B0}=\left(\frac{1}{8}\right)^{2}=\frac{1}{64}
$$

である。

全圧力の値は

$$
(p+p_B)/p_{B0}=\left(1\right) \left(\frac{1}{32}\right) + \left(\frac{1}{64}\right)=\frac{3}{64}
$$

である。

圧縮後のベータは

$$
\beta=\frac{\left(1\right) \left(\frac{1}{32}\right)}{\left(\frac{1}{64}\right)}=2
$$

である。

音速の二乗は

$$
c_s^2/v_{A0}^2=\frac{\left(1\right) \left(\frac{1}{32}\right) \left(\frac{5}{3}\right)}{2 \left(\frac{1}{8}\right)}=\frac{5}{24}
$$

である。

垂直速波の速度の二乗は

$$
c_f^2/v_{A0}^2=\left(\frac{1}{8}\right) + \left(\frac{5}{24}\right)=\frac{1}{3}
$$

である。

この指数では磁気圧が先に下がり、ベータは増える。
:::


## 磁気凍結と磁力線の保存

動く物質曲面 $\Sigma(t)$ に対し、短時間に縁が掃く面の磁束まで足す。Faraday 則と $\nabla\cdot\vb B=0$ から
$$
\dv{}{t}\int_{\Sigma(t)}\vb B\cdot\dd\vb S
=-\oint_{\partial\Sigma(t)}(\vb E+\vb U\times\vb B)\cdot\dd\vb l .
\label{eq-plasma-flux}
$$
右辺がゼロであることを :term[磁気凍結]{en="flux freezing" yomi="じきとうけつ"} という。

### 線要素による証明

誘導方程式を展開して $D_t\vb B=(\vb B\cdot\nabla)\vb U-\vb B\nabla\cdot\vb U$ を得る。連続式から $D_t(\vb B/\rho)=(\vb B/\rho\cdot\nabla)\vb U$ である。流体粒子間の微小線要素 $\delta\vb x$ も $D_t\delta\vb x=(\delta\vb x\cdot\nabla)\vb U$ を満たす。同じ線形常微分方程式の一意性により初めに平行なら平行のままである。滑らかで可逆な流れ写像の間は、磁力線の連結関係が保たれる。

:::example[横圧縮]
$x$ 方向の長さだけが半分になり $\vb B\parallel\hat{\vb z}$ なら、$z$ を法線とする物質断面積も半分になる。従って $B$ と $\rho$ はともに二倍となる。
:::

### 流れ写像と Cauchy の公式

初期ラベル $\vb a$ から位置 $\vb x(\vb a,t)$ への写像の勾配を $G=\partial\vb x/\partial\vb a$ とする。鎖律から $\dot G=(\nabla\vb U)G$、Jacobi の行列式公式から $\dot J=J\nabla\cdot\vb U$、$J=\det G$ である。密度は $\rho=\rho_0/J$ となる。磁場の式に代入すると
$$
\vb B(\vb x(\vb a,t),t)=\frac{G(\vb a,t)\vb B_0(\vb a)}{J(\vb a,t)}.
$$
実際、$\vb B/\rho$ と $G\vb B_0/\rho_0$ は同じ線形初期値問題を満たす。面積ベクトルは $J G^{-\mathsf T}$ で変換されるので、$\vb B$ との内積において $J$ と $G$ が消え、各微小面積の磁束が不変になる。

:::example[せん断で生まれる磁場]
$\vb U=(a y,0,0)$、$\vb B_0=B_0\hat{\vb y}$ なら $\vb x=(a_x+at a_y,a_y,a_z)$ で $J=1$ である。よって $\vb B=B_0(at,1,0)$ となる。磁束が保存されても、同じ空間座標での磁場強度は変化できる。
:::

### せん断する流れ写像と磁力線の傾き

流れ $\vb u=a y\hat{\vb x}$ とし密度一定とする。写像は $(X,Y,Z)\mapsto(X+atY,Y,Z)$ で行列式は一である。変形勾配は $I+at\,\hat{\vb x}\hat{\vb y}$ だから Cauchy の公式より $B_x=B_{x0}+at B_{y0}$ となる。$B_y,B_z$ は一定である。初めに $xy$ 面内にある磁力線に沿う微小変位も全く同じ行列を掛けられるので、変形後の磁場と平行なまま流れとともに移る。磁場のエネルギーが増えるときは、流れとの仕事や境界の Poynting 流束を含めて収支を取る。ここでは流れを与えて誘導式だけを解いている。運動方程式も課す場合は、磁気力と境界条件を別に確認する。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $a$ | せん断率 | $\mathrm{s^{-1}}$ |
| $t$ | 経過時間 | $\mathrm s$ |
| $B_{x0}/B_*$ | 初期の横向き成分 | 1 |
| $B_{y0}/B_*$ | 初期のせん断を横切る成分 | 1 |

せん断の累積値を $s$ とする。流れ写像の非対角成分である。

$$
s=a t
$$

増幅後の横成分を $B_x/B_*$ とする。初期の $B_y$ が横方向へ伸ばされる。

$$
B_x/B_*=B_{x0}/B_* + B_{y0}/B_* s
$$

磁気エネルギーの規格化値を $B^2/B_*^2$ とする。$B_z=0$ を仮定する。

$$
B^2/B_*^2=(B_x/B_*)^2 + (B_{y0}/B_*)^2
$$

初期と比べたエネルギーを $B^2/B_0^2$ とする。圧縮のない変形でも変わる。

$$
B^2/B_0^2=\frac{B^2/B_*^2}{(B_{x0}/B_*)^2 + (B_{y0}/B_*)^2}
$$

磁力線の局所傾きを $B_x/B_y$ とする。せん断方向への傾斜に相当する。

$$
B_x/B_y=\frac{(B_x/B_*)}{B_{y0}/B_*}
$$

瞬時の磁気エネルギー増加率を $\mu_0\dot U_B/B_*^2$ とする。$U_B=B^2/(2\mu_0)$ を微分する。

$$
\mu_0\dot U_B/B_*^2=(B_x/B_*) B_{y0}/B_* a
$$

:::example[直交する初期磁場]
初期に $B_y$ だけがあるときの増幅を求める。

入力値を次のように置く。

$$
\begin{aligned}
a&=1\,\mathrm{s^{-1}} \\t&=2\,\mathrm s \\B_{x0}/B_*&=0 \\B_{y0}/B_*&=1
\end{aligned}
$$

これらを順に代入する。

せん断の累積値は

$$
s=\left(1\right) \left(2\right)=2
$$

である。

増幅後の横成分は

$$
B_x/B_*=\left(0\right) + \left(1\right) \left(2\right)=2
$$

である。

磁気エネルギーの規格化値は

$$
B^2/B_*^2=\left(1\right)^{2} + \left(2\right)^{2}=5
$$

である。

初期と比べたエネルギーは

$$
B^2/B_0^2=\frac{\left(5\right)}{\left(0\right)^{2} + \left(1\right)^{2}}=5
$$

である。

磁力線の局所傾きは

$$
B_x/B_y=\frac{\left(2\right)}{\left(1\right)}=2
$$

である。

瞬時の磁気エネルギー増加率は

$$
\mu_0\dot U_B/B_*^2=\left(1\right)^{2} \left(2\right)=2\,\mathrm{s^{-1}}
$$

である。

磁力線は流れに沿って傾き、エネルギーは二次で増す。
:::

:::example[初期の傾きを戻す時間]
せん断が最初の傾きを打ち消す場合を求める。

入力値を次のように置く。

$$
\begin{aligned}
a&=1\,\mathrm{s^{-1}} \\t&=1\,\mathrm s \\B_{x0}/B_*&=-1 \\B_{y0}/B_*&=1
\end{aligned}
$$

これらを順に代入する。

せん断の累積値は

$$
s=\left(1\right)^{2}=1
$$

である。

増幅後の横成分は

$$
B_x/B_*=\left(-1\right) + \left(1\right)^{2}=0
$$

である。

磁気エネルギーの規格化値は

$$
B^2/B_*^2=\left(0\right)^{2} + \left(1\right)^{2}=1
$$

である。

初期と比べたエネルギーは

$$
B^2/B_0^2=\frac{\left(1\right)}{\left(-1\right)^{2} + \left(1\right)^{2}}=\frac{1}{2}
$$

である。

磁力線の局所傾きは

$$
B_x/B_y=\frac{\left(0\right)}{\left(1\right)}=0
$$

である。

瞬時の磁気エネルギー増加率は

$$
\mu_0\dot U_B/B_*^2=\left(0\right) \left(1\right)^{2}=0\,\mathrm{s^{-1}}
$$

である。

この時刻に磁気エネルギーは極小となる。
:::

:::example[一時的に減衰する磁場]
磁気拡散を入れなくてもエネルギーが減るか調べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=1\,\mathrm{s^{-1}} \\t&=1\,\mathrm s \\B_{x0}/B_*&=-2 \\B_{y0}/B_*&=1
\end{aligned}
$$

これらを順に代入する。

せん断の累積値は

$$
s=\left(1\right)^{2}=1
$$

である。

増幅後の横成分は

$$
B_x/B_*=\left(-2\right) + \left(1\right)^{2}=-1
$$

である。

磁気エネルギーの規格化値は

$$
B^2/B_*^2=\left(-1\right)^{2} + \left(1\right)^{2}=2
$$

である。

初期と比べたエネルギーは

$$
B^2/B_0^2=\frac{\left(2\right)}{\left(-2\right)^{2} + \left(1\right)^{2}}=\frac{2}{5}
$$

である。

磁力線の局所傾きは

$$
B_x/B_y=\frac{\left(-1\right)}{\left(1\right)}=-1
$$

である。

瞬時の磁気エネルギー増加率は

$$
\mu_0\dot U_B/B_*^2=\left(-1\right) \left(1\right)^{2}=-1\,\mathrm{s^{-1}}
$$

である。

せん断が磁場を短くする段階であり、凍結条件と両立する。
:::

:::example[逆向きのせん断]
増幅率と磁力線の傾きの符号を比べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=-1\,\mathrm{s^{-1}} \\t&=2\,\mathrm s \\B_{x0}/B_*&=0 \\B_{y0}/B_*&=1
\end{aligned}
$$

これらを順に代入する。

せん断の累積値は

$$
s=\left(-1\right) \left(2\right)=-2
$$

である。

増幅後の横成分は

$$
B_x/B_*=\left(-2\right) \left(1\right) + \left(0\right)=-2
$$

である。

磁気エネルギーの規格化値は

$$
B^2/B_*^2=\left(-2\right)^{2} + \left(1\right)^{2}=5
$$

である。

初期と比べたエネルギーは

$$
B^2/B_0^2=\frac{\left(5\right)}{\left(0\right)^{2} + \left(1\right)^{2}}=5
$$

である。

磁力線の局所傾きは

$$
B_x/B_y=\frac{\left(-2\right)}{\left(1\right)}=-2
$$

である。

瞬時の磁気エネルギー増加率は

$$
\mu_0\dot U_B/B_*^2=\left(-1\right) \left(-2\right) \left(1\right)=2\,\mathrm{s^{-1}}
$$

である。

初期に $B_y$ だけなら符号反転は傾きの向きを変え、エネルギーは同じである。
:::


## 磁気圧と磁気張力

磁気力は場の強度の勾配と、向きが曲がる効果に分けられる。$\nabla\cdot\vb B=0$ で成分の外積を展開すると
$$
\vb j\times\vb B=-\nabla(B^2/2\mu_0)+(\vb B\cdot\nabla)\vb B/\mu_0
=\nabla\cdot\left(\frac{\vb B\vb B}{\mu_0}-\frac{B^2}{2\mu_0}I\right).
\label{eq-plasma-magnetic-stress}
$$
第 1 項の単位は圧力の勾配である。第 2 項の垂直成分は $B^2\vb\kappa/\mu_0$ であり、磁力線の曲率中心へ向く。平行成分は二項で打ち消し合う。

:::example[円形磁場]
$\vb B=B_\theta(r)\hat{\vb\theta}$ なら力の動径成分は $-\partial_r(B_\theta^2/(2\mu_0))-B_\theta^2/(\mu_0 r)$ である。$B_\theta=C/r$ で二項が相殺し、この領域の電流はゼロである。
:::

:::example[張力としての境界応力]
磁気応力 $\mathsf T_B=\vb B\vb B/\mu_0-I B^2/(2\mu_0)$ を、法線 $\hat{\vb n}\parallel\vb B$ の面に作用させると $+\hat{\vb n}B^2/(2\mu_0)$ である。$\hat{\vb n}\perp\vb B$ なら $-\hat{\vb n}B^2/(2\mu_0)$ である。法線方向を変えることで場の引張と圧力の違いを読める。
:::

### 円筒磁場で磁気圧と張力を分解する

$\vb B=B_\theta(r)\hat{\vb\theta}$、$B_\theta=B_a(r/a)^p$ とする。$\vb B\cdot\grad\vb B$ の動径成分は基底の微分により $-B_\theta^2/r$ である。磁気圧の微分は $-\partial_r(B_\theta^2/2\mu_0)=-p B_\theta^2/(\mu_0r)$ なので、和は $-(p+1)B_\theta^2/(\mu_0r)$ となる。Ampère の法則でも $J_z=(r B_\theta)'/(\mu_0r)$ となり、$-\!J_z B_\theta$ と一致する。この検査は円筒軸から離れた場所で行う。負のべきの磁場を軸まで延ばすと特異になる場合がある。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $x$ | 動径 $r/a>0$ | 1 |
| $p$ | 磁場分布のべき指数 | 1 |

磁場の比を $B_\theta/B_a$ とする。局所での幾何学を決める。

$$
B_\theta/B_a=x^{p}
$$

磁気張力を $F_t/F_*$ とする。力密度の単位は $F_*=B_a^2/(\mu_0a)$ である。

$$
F_t/F_*=- \frac{(B_\theta/B_a)^{2}}{x}
$$

磁気圧による力密度を $F_p/F_*$ とする。強度の勾配の符号で向きが変わる。

$$
F_p/F_*=(F_t/F_*) p
$$

磁気力密度を $F_r/F_*$ とする。張力と圧力のベクトル和である。

$$
F_r/F_*=(F_p/F_*) + (F_t/F_*)
$$

軸電流密度を $\mu_0aJ_z/B_a$ とする。Ampère の法則を円筒形で使う。

$$
\mu_0aJ_z/B_a=\frac{(B_\theta/B_a) \left(p + 1\right)}{x}
$$

静止平衡に必要な圧力勾配を $p'(r)/F_*$ とする。$-\grad p+\vb J\times\vb B=0$ による。

$$
p'(r)/F_*=(F_r/F_*)
$$

:::example[一様軸電流]
半径に比例する磁場から力を計算する。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{2} \\p&=1
\end{aligned}
$$

これらを順に代入する。

磁場の比は

$$
B_\theta/B_a=\left(\frac{1}{2}\right)^{\left(1\right)}=\frac{1}{2}
$$

である。

磁気張力は

$$
F_t/F_*=- \left(\frac{1}{2}\right)=- \frac{1}{2}
$$

である。

磁気圧による力密度は

$$
F_p/F_*=\left(- \frac{1}{2}\right) \left(1\right)=- \frac{1}{2}
$$

である。

磁気力密度は

$$
F_r/F_*=2 \left(- \frac{1}{2}\right)=-1
$$

である。

軸電流密度は

$$
\mu_0aJ_z/B_a=\left(1\right) + 1=2
$$

である。

静止平衡に必要な圧力勾配は

$$
p'(r)/F_*=\left(-1\right)=-1
$$

である。

圧力と張力が同じ大きさで内向きに働く。
:::

:::example[真空外部の磁場]
電流のない周囲の二つの力を求める。

入力値を次のように置く。

$$
\begin{aligned}
x&=2 \\p&=-1
\end{aligned}
$$

これらを順に代入する。

磁場の比は

$$
B_\theta/B_a=\left(2\right)^{\left(-1\right)}=\frac{1}{2}
$$

である。

磁気張力は

$$
F_t/F_*=- \frac{\left(\frac{1}{2}\right)^{2}}{\left(2\right)}=- \frac{1}{8}
$$

である。

磁気圧による力密度は

$$
F_p/F_*=\left(- \frac{1}{8}\right) \left(-1\right)=\frac{1}{8}
$$

である。

磁気力密度は

$$
F_r/F_*=\left(- \frac{1}{8}\right) + \left(\frac{1}{8}\right)=0
$$

である。

軸電流密度は

$$
\mu_0aJ_z/B_a=\frac{\left(\frac{1}{2}\right) \left(\left(-1\right) + 1\right)}{\left(2\right)}=0
$$

である。

静止平衡に必要な圧力勾配は

$$
p'(r)/F_*=\left(0\right)=0
$$

である。

各項は有限だが打ち消し、Lorentz 力は消える。
:::

:::example[一定強度の円形磁場]
強度の勾配をなくした場合を考える。

入力値を次のように置く。

$$
\begin{aligned}
x&=1 \\p&=0
\end{aligned}
$$

これらを順に代入する。

磁場の比は

$$
B_\theta/B_a=\left(1\right)^{\left(0\right)}=1
$$

である。

磁気張力は

$$
F_t/F_*=- \left(1\right)=-1
$$

である。

磁気圧による力密度は

$$
F_p/F_*=\left(-1\right) \left(0\right)=0
$$

である。

磁気力密度は

$$
F_r/F_*=\left(-1\right) + \left(0\right)=-1
$$

である。

軸電流密度は

$$
\mu_0aJ_z/B_a=\left(0\right) + 1=1
$$

である。

静止平衡に必要な圧力勾配は

$$
p'(r)/F_*=\left(-1\right)=-1
$$

である。

張力は残るため、平衡には減少する熱圧力が必要である。
:::

:::example[急減する環状場]
局所的に正の動径力を作る例を調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=2 \\p&=-2
\end{aligned}
$$

これらを順に代入する。

磁場の比は

$$
B_\theta/B_a=\left(2\right)^{\left(-2\right)}=\frac{1}{4}
$$

である。

磁気張力は

$$
F_t/F_*=- \frac{\left(\frac{1}{4}\right)^{2}}{\left(2\right)}=- \frac{1}{32}
$$

である。

磁気圧による力密度は

$$
F_p/F_*=\left(- \frac{1}{32}\right) \left(-2\right)=\frac{1}{16}
$$

である。

磁気力密度は

$$
F_r/F_*=\left(- \frac{1}{32}\right) + \left(\frac{1}{16}\right)=\frac{1}{32}
$$

である。

軸電流密度は

$$
\mu_0aJ_z/B_a=\frac{\left(\frac{1}{4}\right) \left(\left(-2\right) + 1\right)}{\left(2\right)}=- \frac{1}{8}
$$

である。

静止平衡に必要な圧力勾配は

$$
p'(r)/F_*=\left(\frac{1}{32}\right)=\frac{1}{32}
$$

である。

ここでは外向きの磁気圧力が張力を上回る。
:::


## Alfvén 波と磁気音波

一様静止平衡の磁場を $B_0\hat{\vb z}$、波数を $(k_\perp,0,k_\parallel)$ とする。連続式と圧力式から $\delta p=\rho_0 c_s^2\vb k\cdot\vb u/\omega$、$c_s^2=\gamma p_0/\rho_0$、誘導式から $\delta\vb B=[\vb B_0(\vb k\cdot\vb u)-\vb u(\vb k\cdot\vb B_0)]/\omega$ である。運動方程式の $y$ 成分には圧力が無く、
$$
\omega^2=k_\parallel^2 v_A^2,\qquad v_A=B_0/\sqrt{\mu_0\rho_0}
\label{eq-plasma-alfven}
$$
が :term[Alfvén 波]{en="Alfvén wave" yomi="あるふべんは"} の分散である。$xz$ 成分の運動式は
$$
\omega^2\begin{pmatrix}u_x\\u_z\end{pmatrix}
=\begin{pmatrix}
v_A^2k^2+c_s^2 k_\perp^2 & c_s^2 k_\perp k_\parallel\\
c_s^2 k_\perp k_\parallel & c_s^2 k_\parallel^2
\end{pmatrix}\begin{pmatrix}u_x\\u_z\end{pmatrix}.
$$
従って $\omega^4-k^2(v_A^2+c_s^2)\omega^2+k^2 k_\parallel^2 v_A^2 c_s^2=0$ である。正の平方根の速い枝と遅い枝が :term[磁気音波]{en="magnetosonic wave" yomi="じきおんぱ"} である。垂直方向では遅い枝の伝播速度がゼロになり、速い枝は $\sqrt{c_s^2+v_A^2}$ になる。

:::example[平行伝播の縮退]
$k_\perp=0$ なら $xz$ 行列は対角で、位相速度は $v_A$ と $c_s$ である。前者は $y$ 偏波の Alfvén 波と縮退する。
:::

### 一様平衡のエネルギー分配

平行伝播の Alfvén 波で $\omega=k_\parallel v_A>0$ とすると $\delta\vb B=-\sqrt{\mu_0\rho_0}\vb u_\perp$ である。平均の運動エネルギー密度 $\rho_0|\vb u|^2/4$ と磁場のエネルギー密度 $|\delta\vb B|^2/(4\mu_0)$ が等しい。群速度は $\partial_{\vb k}\omega=v_A\vb b$ である。背景の磁場に斜めの波数でも Alfvén 分枝のエネルギーは磁場に沿って運ばれる。

:::example[斜めの位相面]
$k_\perp=3 k_\parallel$、$k_\parallel>0$ なら Alfvén 波の位相速度の大きさは $v_A/\sqrt{10}$ である。群速度は $v_A$ のまま磁場に沿う。位相面に垂直な方向とエネルギーの移動方向を区別する。
:::

### 非圧縮波を反対向きに分ける

密度一定の非圧縮理想 MHD では、磁場を速度と同じ単位で $\vb B/\sqrt{\mu_0\rho_0}=v_A\hat{\vb z}+\vb b_1$ と書ける。運動式では磁気圧を熱圧力と合わせて全圧力 $\Pi=(p+B^2/(2\mu_0))/\rho_0$ に入れる。運動式と誘導式を足し引きし、:term[Elsasser 変数]{en="Elsasser variables" yomi="えるさっさーへんすう"} $\vb z^\pm=\vb U\pm\vb b_1$ で整理すると
$$
\partial_t\vb z^\pm\mp v_A\partial_z\vb z^\pm+
(\vb z^\mp\cdot\nabla)\vb z^\pm=-\nabla\Pi,qquad
\nabla\cdot\vb z^\pm=0 .
$$
片方がゼロなら、全圧力を一様に取れる周期箱で他方は形を保ったまま場に沿って伝わる。反対向きに進む二つの成分がそろうと、移流項が波形を変える。各式と $\vb z^\pm$ の内積を体積で積分すると、移流項と圧力項は周期境界で発散積分として消える。従って $\int|\vb z^\pm|^2\dd V$ がそれぞれ保存される。その差を四で割ると、揺らぎについての交差ヘリシティ $\int\vb U\cdot\vb b_1\dd V$ になる。

### 速波と遅波を波面に沿って選別する

位相速度の二乗 $X=(\omega/kv_A)^2$ は $X^2-(1+c_s^2/v_A^2)X+(c_s^2/v_A^2)\cos^2\theta=0$ を満たす。縦と磁場に垂直な変位の二成分の行列を対角化した形である。判別式は $(1-c_s^2/v_A^2)^2+4(c_s^2/v_A^2)\sin^2\theta$ で非負になる。低い方の根を根の積から出せば、垂直伝播や小音速での差し引きを避けられる。垂直の遅波の根はゼロであり、伝播する有限周波数枝を意味しない。剪断 Alfvén 波の偏波は波数と磁場の作る平面から外れ、速度は $v_A|\cos\theta|$ である。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $A$ | 音速と Alfvén 速度の二乗比 $c_s^2/v_A^2$ | 1 |
| $C$ | $\cos^2\theta$ | 1 |

判別式を $D$ とする。正の二つの項に分けられる。

$$
D=4 A \left(1 - C\right) + \left(1 - A\right)^{2}
$$

速波の二乗速度を $X_f$ とする。大きな根を選ぶ。

$$
X_f=\frac{A}{2} + \frac{\sqrt{D}}{2} + \frac{1}{2}
$$

遅波の二乗速度を $X_s$ とする。二根の積は $AC$ である。

$$
X_s=\frac{A C}{X_{f}}
$$

速波の位相速度比を $v_f/v_A$ とする。波面法線方向の速度である。

$$
v_f/v_A=\sqrt{X_{f}}
$$

遅波の位相速度比を $v_s/v_A$ とする。垂直入射でゼロに近づく。

$$
v_s/v_A=\sqrt{X_{s}}
$$

剪断 Alfvén 波の位相速度比を $v_{\rm shear}/v_A$ とする。磁場に平行な波数が復元に寄与する。

$$
v_{\rm shear}/v_A=\sqrt{C}
$$

:::example[冷たい平行伝播]
音波と Alfvén 波の速度を根から識別する。

入力値を次のように置く。

$$
\begin{aligned}
A&=\frac{1}{4} \\C&=1
\end{aligned}
$$

これらを順に代入する。

判別式は

$$
D=4 \left(\frac{1}{4}\right) \left(1 - \left(1\right)\right) + \left(1 - \left(\frac{1}{4}\right)\right)^{2}=\frac{9}{16}
$$

である。

速波の二乗速度は

$$
X_f=\frac{\left(\frac{1}{4}\right)}{2} + \frac{\sqrt{\left(\frac{9}{16}\right)}}{2} + \frac{1}{2}=1
$$

である。

遅波の二乗速度は

$$
X_s=\left(\frac{1}{4}\right)=\frac{1}{4}
$$

である。

速波の位相速度比は

$$
v_f/v_A=\sqrt{\left(1\right)}=1
$$

である。

遅波の位相速度比は

$$
v_s/v_A=\sqrt{\left(\frac{1}{4}\right)}=\frac{1}{2}
$$

である。

剪断 Alfvén 波の位相速度比は

$$
v_{\rm shear}/v_A=\sqrt{\left(1\right)}=1
$$

である。

この例では速波と剪断波の速度が一致する。
:::

:::example[高い音速の平行伝播]
音速を増やして縮退する枝を求める。

入力値を次のように置く。

$$
\begin{aligned}
A&=4 \\C&=1
\end{aligned}
$$

これらを順に代入する。

判別式は

$$
D=4 \left(4\right) \left(1 - \left(1\right)\right) + \left(1 - \left(4\right)\right)^{2}=9
$$

である。

速波の二乗速度は

$$
X_f=\frac{\left(4\right)}{2} + \frac{\sqrt{\left(9\right)}}{2} + \frac{1}{2}=4
$$

である。

遅波の二乗速度は

$$
X_s=\left(1\right)=1
$$

である。

速波の位相速度比は

$$
v_f/v_A=\sqrt{\left(4\right)}=2
$$

である。

遅波の位相速度比は

$$
v_s/v_A=\sqrt{\left(1\right)}=1
$$

である。

剪断 Alfvén 波の位相速度比は

$$
v_{\rm shear}/v_A=\sqrt{\left(1\right)}=1
$$

である。

遅波の方が剪断波と同じ速度になる。
:::

:::example[垂直方向]
三つの速度を比較する。

入力値を次のように置く。

$$
\begin{aligned}
A&=1 \\C&=0
\end{aligned}
$$

これらを順に代入する。

判別式は

$$
D=4 \left(1\right) \left(1 - \left(0\right)\right) + \left(1 - \left(1\right)\right)^{2}=4
$$

である。

速波の二乗速度は

$$
X_f=\frac{\left(1\right)}{2} + \frac{\sqrt{\left(4\right)}}{2} + \frac{1}{2}=2
$$

である。

遅波の二乗速度は

$$
X_s=\frac{\left(0\right) \left(1\right)}{\left(2\right)}=0
$$

である。

速波の位相速度比は

$$
v_f/v_A=\sqrt{\left(2\right)}\simeq 1.4142
$$

である。

遅波の位相速度比は

$$
v_s/v_A=\sqrt{\left(0\right)}=0
$$

である。

剪断 Alfvén 波の位相速度比は

$$
v_{\rm shear}/v_A=\sqrt{\left(0\right)}=0
$$

である。

速波の圧縮だけが有限の伝播速度を持つ。
:::

:::example[斜めで音速と Alfvén 速さが等しい]
斜め方向での縮退の解除を確かめる。

入力値を次のように置く。

$$
\begin{aligned}
A&=1 \\C&=\frac{3}{4}
\end{aligned}
$$

これらを順に代入する。

判別式は

$$
D=4 \left(1\right) \left(1 - \left(\frac{3}{4}\right)\right) + \left(1 - \left(1\right)\right)^{2}=1
$$

である。

速波の二乗速度は

$$
X_f=\frac{\sqrt{\left(1\right)}}{2} + \frac{\left(1\right)}{2} + \frac{1}{2}=\frac{3}{2}
$$

である。

遅波の二乗速度は

$$
X_s=\frac{\left(1\right) \left(\frac{3}{4}\right)}{\left(\frac{3}{2}\right)}=\frac{1}{2}
$$

である。

速波の位相速度比は

$$
v_f/v_A=\sqrt{\left(\frac{3}{2}\right)}\simeq 1.2247
$$

である。

遅波の位相速度比は

$$
v_s/v_A=\sqrt{\left(\frac{1}{2}\right)}\simeq 0.70711
$$

である。

剪断 Alfvén 波の位相速度比は

$$
v_{\rm shear}/v_A=\sqrt{\left(\frac{3}{4}\right)}\simeq 0.86603
$$

である。

判別式が有限となり、速波と遅波が分かれる。
:::


## 抵抗性 MHD と磁気 Reynolds 数

抵抗率 $\eta$ を残して $\vb E+\vb U\times\vb B=\eta\vb j$ とする。一様 $\eta$ なら Faraday 則に代入し
$$
\partial_t\vb B=\nabla\times(\vb U\times\vb B)+\eta_m\nabla^2\vb B,\qquad
\eta_m=\eta/\mu_0 .
\label{eq-plasma-resistive-induction}
$$
$\eta_m$ の単位は $\mathrm{m^2\,s^{-1}}$ である。:term[磁気 Reynolds 数]{en="magnetic Reynolds number" yomi="じきれいのるずすう"} $R_m=UL/\eta_m$ は移流と拡散の比である。$U=v_A$ の場合は Lundquist 数 $S_L=Lv_A/\eta_m$ と呼ぶ。

:::example[磁場の減衰]
流れを $\vb U=0$ に固定して誘導式だけを解くと、$B_y=B_1\cos kx$ の初期値は $e^{-\eta_m k^2t}$ で減る。磁気エネルギーの時間指数は二倍であり、周期境界での体積積分の減少率は $\int\eta j^2\,\dd V$ である。
:::

### 抵抗のある Alfvén 波

平行な横摂動で一様平衡を考える。流体式は $\partial_t\vb u_\perp=B_0\partial_z\vb B_1/(\mu_0\rho_0)$、誘導式は $\partial_t\vb B_1=B_0\partial_z\vb u_\perp+\eta_m\partial_z^2\vb B_1$ である。時間指数 $\lambda$ を取れば $\lambda(\lambda+\eta_m k^2)+k^2v_A^2=0$ となる。弱い拡散で $\eta_m k\ll v_A$ なら振幅は $e^{-\eta_m k^2t/2}$ で減衰する。速度と磁場がエネルギーを分け合うため、静止した媒質の磁気拡散の減衰率と異なる。

:::example[振動できる範囲]
$\eta_m k=2v_A$ で平方根の判別式が消え、二つの指数は $-\eta_m k^2/2$ に重なる。$\eta_m k$ がさらに大きいと指数は両方とも負の実数となる。
:::

### 抵抗加熱と総エネルギー

抵抗を入れると磁場が失う仕事率は $\vb j\cdot\vb E=\vb U\cdot(\vb j\times\vb B)+\eta j^2$ となる。最初の項は流体の力学的仕事へ移る。二つ目は内部エネルギーの源になる。内部エネルギー $u=p/(\gamma-1)$ は、熱伝導と粘性を省けば
$$
\partial_t u+\nabla\cdot(u\vb U)=-p\nabla\cdot\vb U+\eta j^2
$$
を満たす。連続式で密度の時間微分を消せば $D_t(p\rho^{-\gamma})=(\gamma-1)\rho^{-\gamma}\eta j^2$ となる。磁気・運動・内部の三式を足すと、総エネルギーは理想 MHD と同じ形の流束で保存する。ただし Poynting 流束には抵抗を含む実際の電場を使う。空間的に変化する $\eta$ でも仕事率は $\eta j^2$ であり、誘導式の拡散項は $-\nabla\times(\eta\vb j)$ の形で残る。

### 拡散と移流の比較を Fourier モードで行う

一定の $\eta_m=\eta/\mu_0$ と一様な流速 $U$ で、横磁場成分は $\partial_tB+U\partial_xB=\eta_m\partial_x^2B$ を満たす。$Be^{-ikx}$ の空間積分で Fourier 振幅を取ると $\dot B_k=(-ikU-\eta_m k^2)B_k$ である。時間 $L/U$ で測った拡散指数は $(kL)^2/R_m$ になる。したがって装置の $R_m=UL/\eta_m$ が大きくても、高い波数を保持した磁場が凍結するとは限らない。一波数のエネルギーは振幅の二乗で減り、その差は周期領域なら抵抗加熱に変わる。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $R_m$ | 長さ $L$ による磁気 Reynolds 数 | 1 |
| $K$ | 波数比 $kL$ | 1 |
| $s$ | 時間 $Ut/L$ | 1 |

移流で進む位相を $\theta$ とする。振幅を変えずに位相だけを変える。

$$
\theta=K s
$$

抵抗性減衰の指数を $A$ とする。波数の二乗で増える。

$$
A=\frac{K^{2} s}{R_{m}}
$$

磁場の振幅比を $\abs{B_k(t)/B_k(0)}$ とする。位相因子の絶対値は一である。

$$
\abs{B_k(t)/B_k(0)}=e^{- A}
$$

磁気エネルギーの残存率を $U_B(t)/U_B(0)$ とする。二乗のため減衰指数が二倍になる。

$$
U_B(t)/U_B(0)=\abs{B_k(t)/B_k(0)}^2
$$

磁気エネルギーの損失率を $1-U_B(t)/U_B(0)$ とする。境界流束がない場合の加熱である。

$$
1-U_B(t)/U_B(0)=1 - U_B(t)/U_B(0)
$$

一波長尺度での拡散と移流の率の比を $\eta_m k/U$ とする。Fourier 時間率 $\eta_m k^2$ と $Uk$ を比べる。

$$
\eta_m k/U=\frac{K}{R_{m}}
$$

:::example[装置尺度の長波長]
大きな $R_m$ で磁場がほぼ運ばれるか計算する。

入力値を次のように置く。

$$
\begin{aligned}
R_m&=1000 \\K&=1 \\s&=1
\end{aligned}
$$

これらを順に代入する。

移流で進む位相は

$$
\theta=\left(1\right)^{2}=1
$$

である。

抵抗性減衰の指数は

$$
A=\frac{\left(1\right)^{3}}{\left(1000\right)}=\frac{1}{1000}
$$

である。

磁場の振幅比は

$$
\abs{B_k(t)/B_k(0)}=e^{- \left(\frac{1}{1000}\right)}\simeq 0.999
$$

である。

磁気エネルギーの残存率は

$$
U_B(t)/U_B(0)=\left(\left(e^{- \frac{1}{1000}}\right)\right)^{2}\simeq 0.998
$$

である。

磁気エネルギーの損失率は

$$
1-U_B(t)/U_B(0)=1 - \left(e^{- \frac{1}{500}}\right)\simeq 0.001998
$$

である。

一波長尺度での拡散と移流の率の比は

$$
\eta_m k/U=\frac{\left(1\right)}{\left(1000\right)}=\frac{1}{1000}
$$

である。

一移流時間で失う振幅は千分の一程度である。
:::

:::example[同じ装置の細かな構造]
高波数に替えて抵抗性を評価する。

入力値を次のように置く。

$$
\begin{aligned}
R_m&=1000 \\K&=100 \\s&=1
\end{aligned}
$$

これらを順に代入する。

移流で進む位相は

$$
\theta=\left(100\right) \left(1\right)=100
$$

である。

抵抗性減衰の指数は

$$
A=\frac{\left(100\right)^{2} \left(1\right)}{\left(1000\right)}=10
$$

である。

磁場の振幅比は

$$
\abs{B_k(t)/B_k(0)}=e^{- \left(10\right)}\simeq 4.54\times 10^{-5}
$$

である。

磁気エネルギーの残存率は

$$
U_B(t)/U_B(0)=\left(\left(e^{-10}\right)\right)^{2}\simeq 2.0612\times 10^{-9}
$$

である。

磁気エネルギーの損失率は

$$
1-U_B(t)/U_B(0)=1-e^{-20}\simeq 0.999999998
$$

である。

一波長尺度での拡散と移流の率の比は

$$
\eta_m k/U=\frac{\left(100\right)}{\left(1000\right)}=\frac{1}{10}
$$

である。

大域的な $R_m$ が同じでも細部は早く消える。
:::

:::example[拡散と移流が同程度]
一ラジアン分の移流で失われる磁場を求める。

入力値を次のように置く。

$$
\begin{aligned}
R_m&=10 \\K&=10 \\s&=\frac{1}{10}
\end{aligned}
$$

これらを順に代入する。

移流で進む位相は

$$
\theta=\left(10\right) \left(\frac{1}{10}\right)=1
$$

である。

抵抗性減衰の指数は

$$
A=\left(10\right) \left(\frac{1}{10}\right)=1
$$

である。

磁場の振幅比は

$$
\abs{B_k(t)/B_k(0)}=e^{- \left(1\right)}\simeq 0.36788
$$

である。

磁気エネルギーの残存率は

$$
U_B(t)/U_B(0)=\left(\left(e^{-1}\right)\right)^{2}\simeq 0.13534
$$

である。

磁気エネルギーの損失率は

$$
1-U_B(t)/U_B(0)=1 - \left(e^{-2}\right)\simeq 0.86466
$$

である。

一波長尺度での拡散と移流の率の比は

$$
\eta_m k/U=1
$$

である。

モードの移流率と拡散率が等しい。
:::

:::example[長時間の観測]
弱い拡散が蓄積する時間を求める。

入力値を次のように置く。

$$
\begin{aligned}
R_m&=1000 \\K&=1 \\s&=1000
\end{aligned}
$$

これらを順に代入する。

移流で進む位相は

$$
\theta=\left(1000\right) \left(1\right)=1000
$$

である。

抵抗性減衰の指数は

$$
A=\left(1\right)^{2}=1
$$

である。

磁場の振幅比は

$$
\abs{B_k(t)/B_k(0)}=e^{- \left(1\right)}\simeq 0.36788
$$

である。

磁気エネルギーの残存率は

$$
U_B(t)/U_B(0)=\left(\left(e^{-1}\right)\right)^{2}\simeq 0.13534
$$

である。

磁気エネルギーの損失率は

$$
1-U_B(t)/U_B(0)=1 - \left(e^{-2}\right)\simeq 0.86466
$$

である。

一波長尺度での拡散と移流の率の比は

$$
\eta_m k/U=\frac{\left(1\right)}{\left(1000\right)}=\frac{1}{1000}
$$

である。

一移流時間で小さかった効果が有限の損失をもたらす。
:::


## MHD のエネルギー原理

流れのない平衡、固定完全導体の境界で、許される変位 $\vb\xi$ を用いる。境界を固定して表面変位をゼロとし、重力を省く。線形摂動は $\vb Q=\delta\vb B=\nabla\times(\vb\xi\times\vb B)$、$\delta p=-\gamma p\nabla\cdot\vb\xi-\vb\xi\cdot\nabla p$ である。運動式は $\rho\ddot{\vb\xi}=F[\vb\xi]=-\nabla\delta p+(\nabla\times\vb Q)\times\vb B/\mu_0+\vb j\times\vb Q$ となる。実変位で
$$
\delta W=-\frac12\int\vb\xi\cdot F[\vb\xi]\,\dd V
=\frac12\int\left[\frac{Q^2}{\mu_0}+\gamma p(\nabla\cdot\vb\xi)^2+
(\vb\xi\cdot\nabla p)\nabla\cdot\vb\xi+\vb j\cdot(\vb\xi\times\vb Q)\right]\dd V .
\label{eq-plasma-energy-principle}
$$
圧力項は一度発散を部分積分する。磁場の微分を含む項では $(\nabla\times\vb Q)\cdot(\vb\xi\times\vb B)$ を部分積分すると $Q^2$ が得られる。

### 変分の意味

理想 MHD の物質座標で、位置写像 $\vb x=\vb x(\vb a)$ の変形勾配を $G$、Jacobian を $J$ とする。密度は $\rho_0/J$、圧力は $p_0/J^\gamma$、磁場は $G\vb B_0/J$ である。ポテンシャルエネルギーは
$$
V[\vb x]=\int \left(\frac{p_0 J^{1-\gamma}}{\gamma-1}+
\frac{\abs{G\vb B_0}^2}{2\mu_0 J}\right)\dd^3a.
$$
$\delta J=J\nabla\cdot\vb\xi$、$\delta(G\vb B_0)=(G\vb B_0\cdot\nabla)\vb\xi$ を使えば、第一変分は $-\int\vb\xi\cdot(-\nabla p+\vb j\times\vb B)\dd V$ である。平衡ではこれがゼロで、第二混合変分は $-\int\vb\zeta\cdot F[\vb\xi]\dd V$ である。混合微分の交換により $F$ は対称となる。

有限次元の変位展開では、正定値の質量行列で正規化して実対称剛性行列を対角化できる。したがって $\delta W<0$ の試験変位は負の $\omega^2$ をもつモードを保証する。連続系では自己共役作用素の定義域と連続スペクトルを扱う必要があるので、ここでは変分式と有限次元判定を示し、一般のスペクトルに関する定式化は :cite[bellan-plasma] を参照する。

:::example[一様平衡での符号]
$\vb j=0$、$\nabla p=0$ なら二つの平方項だけが残る。従って負の $\delta W$ は得られない。磁場を曲げず体積も変えない変位は零モードになりうる。
:::

:::example[試験変位の Rayleigh 商]
$\vb\xi=\xi_0\hat{\vb x}\cos kz$ を一様磁場 $B_0\hat{\vb z}$ に入れる。周期領域で $\delta W/V=k^2B_0^2\xi_0^2/(4\mu_0)$、質量積分は $\int\rho \xi^2\dd V=V\rho\xi_0^2/2$ である。比 $2\delta W/(\int\rho\xi^2\dd V)=k^2v_A^2$ が Alfvén 分枝の $\omega^2$ を再現する。この例では試験変位自体が固有関数である。
:::

### 試験変位の磁気エネルギーを直接計算する

周期箱内の一様平衡で、実変位を複素振幅 $\hat{\vb\xi}$ の Fourier モードの実部とする。$\delta\vb B=i[(\vb k\cdot\vb B)\hat{\vb\xi}-(\vb k\cdot\hat{\vb\xi})\vb B]$、$\delta p=-i\gamma p\,\vb k\cdot\hat{\vb\xi}$ である。体積平均の $\delta W$ は $|\delta\vb B|^2/(4\mu_0)+\gamma p|\vb k\cdot\hat{\vb\xi}|^2/4$ で非負になる。係数の四は実波の二乗平均とエネルギーの二分の一の積である。$\vb B=B_0\hat{\vb z}$、$\vb k=(k_x,0,k_z)$ では磁気項は $B_0^2[(k_x^2+k_z^2)|\hat\xi_x|^2+k_z^2|\hat\xi_y|^2]$ に比例する。$\xi_z$ は磁場そのものを動かさず、圧縮を通じて寄与する。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $k_x L$ | 横の規格化波数 | 1 |
| $k_z L$ | 平行の規格化波数 | 1 |
| $\xi_x/L$ | 横面内変位 | 1 |
| $\xi_y/L$ | 面外変位 | 1 |
| $\xi_z/L$ | 平行変位 | 1 |
| $C_s$ | 音速と Alfvén 速さの二乗比 | 1 |

規格化圧縮振幅を $G$ とする。体積変化を決める内積である。

$$
G=(\xi_x/L) (k_x L) + (\xi_z/L) (k_z L)
$$

規格化磁気エネルギーを $W_B$ とする。$B_0^2/(4\mu_0)$ で割る。

$$
W_B=(\xi_x/L)^{2} \left((k_x L)^{2} + (k_z L)^{2}\right) + (\xi_y/L)^{2} (k_z L)^{2}
$$

規格化圧縮エネルギーを $W_p$ とする。同じ単位なので係数に音速の二乗比が入る。

$$
W_p=C_{s} G^{2}
$$

規格化した全変分を $W_B+W_p$ とする。一様平衡のため駆動の負項がない。

$$
W_B+W_p=W_{B} + W_{p}
$$

変位の質量規格化ノルムを $M$ とする。共通の密度と実波係数を除く。

$$
M=(\xi_x/L)^{2} + (\xi_y/L)^{2} + (\xi_z/L)^{2}
$$

規格化 Rayleigh 商を $\mathcal R$ とする。単位は $v_A^2/L^2$ である。

$$
\mathcal R=\frac{W_B+W_p}{M}
$$

:::example[面外の剪断]
圧縮のない変位が磁場を曲げるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
k_x L&=1 \\k_z L&=2 \\\xi_x/L&=0 \\\xi_y/L&=1 \\\xi_z/L&=0 \\C_s&=1
\end{aligned}
$$

これらを順に代入する。

規格化圧縮振幅は

$$
G=\left(0\right) \left(1\right) + \left(0\right) \left(2\right)=0
$$

である。

規格化磁気エネルギーは

$$
W_B=\left(0\right)^{2} \left(\left(1\right)^{2} + \left(2\right)^{2}\right) + \left(1\right)^{2} \left(2\right)^{2}=4
$$

である。

規格化圧縮エネルギーは

$$
W_p=\left(0\right)^{2} \left(1\right)=0
$$

である。

規格化した全変分は

$$
W_B+W_p=\left(0\right) + \left(4\right)=4
$$

である。

変位の質量規格化ノルムは

$$
M=2 \left(0\right)^{2} + \left(1\right)^{2}=1
$$

である。

規格化 Rayleigh 商は

$$
\mathcal R=\frac{\left(4\right)}{\left(1\right)}=4
$$

である。

Rayleigh 商は平行波数の二乗になる。
:::

:::example[磁場に沿った圧縮]
圧縮だけの変位に対する変分を求める。

入力値を次のように置く。

$$
\begin{aligned}
k_x L&=1 \\k_z L&=2 \\\xi_x/L&=0 \\\xi_y/L&=0 \\\xi_z/L&=1 \\C_s&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

規格化圧縮振幅は

$$
G=\left(0\right) \left(1\right) + \left(1\right) \left(2\right)=2
$$

である。

規格化磁気エネルギーは

$$
W_B=\left(0\right)^{2} \left(2\right)^{2} + \left(0\right)^{2} \left(\left(1\right)^{2} + \left(2\right)^{2}\right)=0
$$

である。

規格化圧縮エネルギーは

$$
W_p=\left(2\right)^{2} \left(\frac{1}{2}\right)=2
$$

である。

規格化した全変分は

$$
W_B+W_p=\left(0\right) + \left(2\right)=2
$$

である。

変位の質量規格化ノルムは

$$
M=2 \left(0\right)^{2} + \left(1\right)^{2}=1
$$

である。

規格化 Rayleigh 商は

$$
\mathcal R=\frac{\left(2\right)}{\left(1\right)}=2
$$

である。

斜め伝播では試験変位が固有偏波とは限らない。
:::

:::example[平行波数のない横圧縮]
磁気圧と熱圧の変分を比較する。

入力値を次のように置く。

$$
\begin{aligned}
k_x L&=2 \\k_z L&=0 \\\xi_x/L&=1 \\\xi_y/L&=0 \\\xi_z/L&=0 \\C_s&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

規格化圧縮振幅は

$$
G=\left(0\right)^{2} + \left(1\right) \left(2\right)=2
$$

である。

規格化磁気エネルギーは

$$
W_B=\left(0\right)^{4} + \left(1\right)^{2} \left(\left(0\right)^{2} + \left(2\right)^{2}\right)=4
$$

である。

規格化圧縮エネルギーは

$$
W_p=\left(2\right)^{2} \left(\frac{1}{2}\right)=2
$$

である。

規格化した全変分は

$$
W_B+W_p=\left(2\right) + \left(4\right)=6
$$

である。

変位の質量規格化ノルムは

$$
M=2 \left(0\right)^{2} + \left(1\right)^{2}=1
$$

である。

規格化 Rayleigh 商は

$$
\mathcal R=\frac{\left(6\right)}{\left(1\right)}=6
$$

である。

平行曲げがなくても横圧縮で磁気エネルギーが増える。
:::

:::example[両方の変分が消える方向]
波面に沿って磁場方向へ動かす。

入力値を次のように置く。

$$
\begin{aligned}
k_x L&=2 \\k_z L&=0 \\\xi_x/L&=0 \\\xi_y/L&=0 \\\xi_z/L&=1 \\C_s&=1
\end{aligned}
$$

これらを順に代入する。

規格化圧縮振幅は

$$
G=\left(0\right) \left(1\right) + \left(0\right) \left(2\right)=0
$$

である。

規格化磁気エネルギーは

$$
W_B=\left(0\right)^{4} + \left(0\right)^{2} \left(\left(0\right)^{2} + \left(2\right)^{2}\right)=0
$$

である。

規格化圧縮エネルギーは

$$
W_p=\left(0\right)^{2} \left(1\right)=0
$$

である。

規格化した全変分は

$$
W_B+W_p=2 \left(0\right)=0
$$

である。

変位の質量規格化ノルムは

$$
M=2 \left(0\right)^{2} + \left(1\right)^{2}=1
$$

である。

規格化 Rayleigh 商は

$$
\mathcal R=\frac{\left(0\right)}{\left(1\right)}=0
$$

である。

圧縮も磁力線の変形もないので中立である。
:::


### 自己共役な力作用素と成長の意味

実変位の変分式を複素変位へ拡張するには、試験変位を複素共役にして、実の双線形形式を半双線形形式に拡張する。$\langle \xi,\zeta\rangle_\rho=\int\rho\,\vb\xi\cdot\conj{\vb\zeta}\,\dd V$ の質量内積のもとで $A=-\rho^{-1}F$ を定義すれば、$\langle A\xi,\xi\rangle_\rho=2\delta W[\xi]$ である。以下は境界条件と定義域を選んで $A$ が自己共役かつ下に有界になる場合の条件付きの結論である。保存する二次エネルギーは $\tfrac12\|\partial_t\xi\|_\rho^2+\delta W[\xi]$ である。実際、時間微分を取ると $\langle\partial_t^2\xi+A\xi,\partial_t\xi\rangle_\rho$ の実部になり、線形運動式でゼロである。

試験変位の変分が負なら、スペクトル測度に対し $\int\lambda\,\dd\mu_\xi(\lambda)<0$ である。従ってある $\epsilon>0$ の区間 $(-\infty,-\epsilon]$ へのスペクトル射影はゼロではない。初期変位をその射影、初期速度をゼロに選ぶと、この区間での解は $\cosh(t\sqrt{-\lambda})$ を各成分に掛けたものになる。そのノルムは少なくとも $\cosh(t\sqrt\epsilon)$ 倍となり、不安定性が分かる。この議論では負のスペクトルが離散固有値である必要はない。逆に $A\ge0$ なら各スペクトル成分は余弦・正弦で時間発展する。核に属する初期速度を与えると、変位のその成分は時間に比例するため、非負の変分だけから変位ノルムの一様有界性まで結論しない。零から正に離れた下界があれば、保存する二次エネルギーから変位ノルムも抑えられる。自己共役の定義域を構成する問題は前節の文献への案内を参照する。


## ダイナモ理論の概観

抵抗があっても速度の伸張が磁場を保つことがある。与えた非圧縮流で式 :eqref[eq-plasma-resistive-induction] を解く課題を :term[運動学的ダイナモ]{en="kinematic dynamo" yomi="うんどうがくてきだいなも"} という。周期境界では磁気エネルギーの変化は
$$
\dv{}{t}\int \frac{B^2}{2\mu_0}\dd V
=\int\frac{B_iB_j\partial_j U_i}{\mu_0}\dd V
-\frac{\eta_m}{\mu_0}\int \abs{\nabla\times\vb B}^2\dd V .
$$
移流項は体積積分で消え、拡散項は部分積分した。伸張が拡散に勝つ配置が必要だが、$R_m$ の大小だけで成長率は決まらない。磁場が増えると磁気力が流れに作用し、運動学的模型から外れる。以下では流れの対称性で成長が制限される例と、乱流の効果を閉包で表す例を扱う。

:::example[局所の線形伸張]
$\vb U=(ax,-ay,0)$ の理想模型で一様な $B_x$ を置くと $D_t B_x=aB_x$ である。$B_x\propto e^{at}$ は局所での増幅を表す。無限空間でこの流れを置いた例から有限領域での自己励起を直接には結論できない。
:::

### 磁気ヘリシティ

完全導体の固定境界で $\vb B\cdot\hat{\vb n}=0$ とする。単連結な領域で大域的なベクトルポテンシャル $\vb B=\nabla\times\vb A$ をとれるとする。このとき $H_m=\int\vb A\cdot\vb B\dd V$ を :term[磁気ヘリシティ]{en="magnetic helicity" yomi="じきへりしてぃ"} と呼ぶ。$\vb A\mapsto\vb A+\nabla\chi$ による差は $\int\nabla\cdot(\chi\vb B)\dd V=0$ である。$\partial_t\vb A=-\vb E-\nabla\Phi$ と $\partial_t\vb B=-\nabla\times\vb E$ を使い、表面項を落とせるとき $\dot H_m=-2\int\vb E\cdot\vb B\dd V$ である。理想場ではゼロであり、抵抗性 MHD では $-2\int\eta\vb j\cdot\vb B\dd V$ になる。エネルギーの散逸は $\eta j^2$、ヘリシティの変化は $\eta\vb j\cdot\vb B$ であるため、両者は尺度に異なる重みを持つ。

:::example[一つの curl 固有モード]
周期領域で $\nabla\times\vb B=k\vb B$、$k\ne0$ なら $\vb A=\vb B/k$ を選べる。流れをゼロに保つ抵抗模型で磁場は $e^{-\eta_m k^2t}$ で減少する。ヘリシティとエネルギーの比 $H_m/\mathcal E_B=2\mu_0/k$ は一定である。この比の絶対値は $|k|$ が大きいモードほど小さい。
:::

### 二次元の反ダイナモの一例

周期的な二次元領域で、$\partial_z=0$、$\vb U$ と $\vb B$ は面内にあり、$\nabla\cdot\vb U=0$、面内磁場の空間平均はゼロとする。このとき $\vb B=\nabla A\times\hat{\vb z}$ と一価の周期関数 $A$ で表せる。実際、Fourier 成分ごとに $\vb k\cdot\vb B_{\vb k}=0$ であり、$\vb k\ne0$ のみなので各成分は $i\vb k\times\hat{\vb z}$ に比例する。誘導式の面内成分をこの形で書けば
$$
\partial_t A+\vb U\cdot\nabla A=\eta_m\nabla^2 A
$$
を得る。空間一定の項は $A$ の平均をゼロに選べば消える。$A$ を掛けて積分すると $d\int A^2/2\dd V/dt=-\eta_m\int\abs{\nabla A}^2\dd V$ である。最小の非零波数を $k_{\min}$ として Fourier モードごとに足せば $\int\abs{\nabla A}^2\dd V\ge k_{\min}^2\int A^2\dd V$ である。従って $\|A(t)\|_2\le e^{-\eta_m k_{\min}^2t}\|A(0)\|_2$ となる。固定した波形が正の指数で成長する固有モードは存在しない。勾配を増やす流れでは磁場が一時的に増幅することがあり、この不等式は任意の時刻で磁場エネルギーが減るという主張ではない。

:::example[一時的な増幅の解釈]
大域的な磁束ポテンシャルが小さくても、細かい磁場の構造は $A$ の勾配を大きくしうる。例えば $A_j=\cos(jx)/j$ では $\|A_j\|_2$ は $j^{-1}$ に比例するが $\|\partial_x A_j\|_2$ は一定である。反ダイナモの例で排除したのは、一つの固有分布が時間とともに成長を続ける解である。
:::

### 平均場ダイナモの最小模型

平均誘導式に一定の $\alpha$ と全磁気拡散 $\eta_T$ を入れ、$\partial_t\vb B=\alpha\rot\vb B+\eta_T\nabla^2\vb B$ とする。これは乱流起電力を $\alpha\vb B-\eta_t\rot\vb B$ と近似した閉包であり、$\alpha$ が一定であることをここで仮定する。周期箱の螺旋偏波は $\rot\vb B_\pm=\pm k\vb B_\pm$ かつ $\nabla^2\vb B_\pm=-k^2\vb B_\pm$ である。ゆえに $\gamma_\pm=\pm\alpha k-\eta_Tk^2$ となる。成長には一つの偏波で $\alpha k$ が拡散を上回らなければならない。連続の波数における最大点 $k_*=\abs{\alpha}/(2\eta_T)$ が箱に許されるとは限らない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $\alpha$ | 平均起電力の係数 | $\mathrm{m\,s^{-1}}$ |
| $\eta_T$ | 全磁気拡散 | $\mathrm{m^2\,s^{-1}}$ |
| $k$ | 周期箱で許された波数 | $\mathrm{m^{-1}}$ |
| $t$ | 経過時間 | $\mathrm s$ |

正螺旋偏波の成長率を $\gamma_+$ とする。curl 固有値が正である。

$$
\gamma_+=\alpha k - \eta_{T} k^{2}
$$

負螺旋偏波の成長率を $\gamma_-$ とする。起電力だけ符号が逆である。

$$
\gamma_-=- \alpha k - \eta_{T} k^{2}
$$

連続近似の最速波数を $k_{\max}$ とする。二次の成長率を微分する。

$$
k_{\max}=\frac{\left|{\alpha}\right|}{2 \eta_{T}}
$$

連続近似での最大成長率を $\gamma_{\max}$ とする。最速波数を成長率へ戻す。

$$
\gamma_{\max}=\frac{\alpha^{2}}{4 \eta_{T}}
$$

正偏波の振幅増幅率を $G_+$ とする。線形の時間範囲で使う。

$$
G_+=e^{\gamma_{+} t}
$$

負偏波の振幅増幅率を $G_-$ とする。飽和や流れへの反作用を含まない。

$$
G_-=e^{\gamma_{-} t}
$$

:::example[最速波数を選ぶ場合]
二偏波の振幅がどれほど離れるか求める。

入力値を次のように置く。

$$
\begin{aligned}
\alpha&=4\,\mathrm{m\,s^{-1}} \\\eta_T&=1\,\mathrm{m^2\,s^{-1}} \\k&=2\,\mathrm{m^{-1}} \\t&=1\,\mathrm s
\end{aligned}
$$

これらを順に代入する。

正螺旋偏波の成長率は

$$
\gamma_+=- \left(1\right) \left(2\right)^{2} + \left(2\right) \left(4\right)=4\,\mathrm{s^{-1}}
$$

である。

負螺旋偏波の成長率は

$$
\gamma_-=- \left(1\right) \left(2\right)^{2} - \left(2\right) \left(4\right)=-12\,\mathrm{s^{-1}}
$$

である。

連続近似の最速波数は

$$
k_{\max}=\frac{\left|{\left(4\right)}\right|}{2 \left(1\right)}=2\,\mathrm{m^{-1}}
$$

である。

連続近似での最大成長率は

$$
\gamma_{\max}=\frac{\left(4\right)^{2}}{4 \left(1\right)}=4\,\mathrm{s^{-1}}
$$

である。

正偏波の振幅増幅率は

$$
G_+=e^{\left(1\right) \left(4\right)}\simeq 54.598
$$

である。

負偏波の振幅増幅率は

$$
G_-=e^{\left(-12\right) \left(1\right)}\simeq 6.1442\times 10^{-6}
$$

である。

正偏波は成長し、負偏波は拡散と起電力の両方で減る。
:::

:::example[短すぎる波長]
強い起電力でも拡散で減衰する場合を調べる。

入力値を次のように置く。

$$
\begin{aligned}
\alpha&=4\,\mathrm{m\,s^{-1}} \\\eta_T&=1\,\mathrm{m^2\,s^{-1}} \\k&=5\,\mathrm{m^{-1}} \\t&=1\,\mathrm s
\end{aligned}
$$

これらを順に代入する。

正螺旋偏波の成長率は

$$
\gamma_+=- \left(1\right) \left(5\right)^{2} + \left(4\right) \left(5\right)=-5\,\mathrm{s^{-1}}
$$

である。

負螺旋偏波の成長率は

$$
\gamma_-=- \left(1\right) \left(5\right)^{2} - \left(4\right) \left(5\right)=-45\,\mathrm{s^{-1}}
$$

である。

連続近似の最速波数は

$$
k_{\max}=\frac{\left|{\left(4\right)}\right|}{2 \left(1\right)}=2\,\mathrm{m^{-1}}
$$

である。

連続近似での最大成長率は

$$
\gamma_{\max}=\frac{\left(4\right)^{2}}{4 \left(1\right)}=4\,\mathrm{s^{-1}}
$$

である。

正偏波の振幅増幅率は

$$
G_+=e^{\left(-5\right) \left(1\right)}\simeq 0.0067379
$$

である。

負偏波の振幅増幅率は

$$
G_-=e^{\left(-45\right) \left(1\right)}\simeq 2.8625\times 10^{-20}
$$

である。

拡散が波数の二乗で増えるので両偏波が減衰する。
:::

:::example[起電力の符号反転]
螺旋偏波の役割が入れ替わるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
\alpha&=-4\,\mathrm{m\,s^{-1}} \\\eta_T&=1\,\mathrm{m^2\,s^{-1}} \\k&=2\,\mathrm{m^{-1}} \\t&=1\,\mathrm s
\end{aligned}
$$

これらを順に代入する。

正螺旋偏波の成長率は

$$
\gamma_+=\left(-4\right) \left(2\right) - \left(1\right) \left(2\right)^{2}=-12\,\mathrm{s^{-1}}
$$

である。

負螺旋偏波の成長率は

$$
\gamma_-=- \left(-4\right) \left(2\right) - \left(1\right) \left(2\right)^{2}=4\,\mathrm{s^{-1}}
$$

である。

連続近似の最速波数は

$$
k_{\max}=\frac{\left|{\left(-4\right)}\right|}{2 \left(1\right)}=2\,\mathrm{m^{-1}}
$$

である。

連続近似での最大成長率は

$$
\gamma_{\max}=\frac{\left(-4\right)^{2}}{4 \left(1\right)}=4\,\mathrm{s^{-1}}
$$

である。

正偏波の振幅増幅率は

$$
G_+=e^{\left(-12\right) \left(1\right)}\simeq 6.1442\times 10^{-6}
$$

である。

負偏波の振幅増幅率は

$$
G_-=e^{\left(1\right) \left(4\right)}\simeq 54.598
$$

である。

最大成長率は同じで、成長する偏波だけが反転する。
:::

:::example[小さい箱の境界例]
箱の最小波数がちょうど成長限界のときを考える。

入力値を次のように置く。

$$
\begin{aligned}
\alpha&=2\,\mathrm{m\,s^{-1}} \\\eta_T&=1\,\mathrm{m^2\,s^{-1}} \\k&=2\,\mathrm{m^{-1}} \\t&=1\,\mathrm s
\end{aligned}
$$

これらを順に代入する。

正螺旋偏波の成長率は

$$
\gamma_+=- \left(1\right) \left(2\right)^{2} + \left(2\right)^{2}=0\,\mathrm{s^{-1}}
$$

である。

負螺旋偏波の成長率は

$$
\gamma_-=- \left(1\right) \left(2\right)^{2} - \left(2\right)^{2}=-8\,\mathrm{s^{-1}}
$$

である。

連続近似の最速波数は

$$
k_{\max}=\frac{\left|{\left(2\right)}\right|}{2 \left(1\right)}=1\,\mathrm{m^{-1}}
$$

である。

連続近似での最大成長率は

$$
\gamma_{\max}=\frac{\left(2\right)^{2}}{4 \left(1\right)}=1\,\mathrm{s^{-1}}
$$

である。

正偏波の振幅増幅率は

$$
G_+=e^{\left(0\right) \left(1\right)}=1
$$

である。

負偏波の振幅増幅率は

$$
G_-=e^{\left(-8\right) \left(1\right)}\simeq 0.00033546
$$

である。

連続波数で最適な値が存在しても箱内では正成長を得られない。
:::


## 演習問題

::::exercise{level=basic}
$B$ を二倍、密度を四倍にしたとき $v_A$ を比較せよ。

:::solution
$B/\sqrt{\mu_0\rho}$ なので倍率は $2/\sqrt4=1$ である。
:::
::::

::::exercise{level=basic}
$L=2\,\mathrm m,U=3\,\mathrm{m\,s^{-1}},\eta_m=0.1\,\mathrm{m^2\,s^{-1}}$ の $R_m$ を求めよ。

:::solution
$UL/\eta_m=6/0.1=60$ である。
:::
::::

::::exercise{level=basic}
垂直速波で $v_A=3$、$c_s=4$ を同じ速度単位で与える。位相速度は何か。

:::solution
二乗の和の平方根 $\sqrt{9+16}=5$ である。
:::
::::

::::exercise{level=basic}
速度ゼロで $\eta_m=2$、$k=3$ を SI で与える。磁場振幅の減衰率を求めよ。

:::solution
$\eta_m k^2=2\cdot9=18\,\mathrm{s^{-1}}$ である。
:::
::::

::::exercise{level=advanced}
理想 MHD の流速を $\vb U=(ax,by,cz)$ と与えたとき、流体粒子に沿う $B_z$ の変化を求めよ。

:::solution
$\nabla\cdot\vb U=a+b+c$、$[(\vb B\cdot\nabla)\vb U]_z=cB_z$ なので $D_tB_z=-(a+b)B_z$ である。よって $B_z=B_z(0)e^{-(a+b)t}$ となる。
:::
::::

::::exercise{level=advanced}
$B_\theta=\alpha r$ の磁気力を求めよ。

:::solution
磁気圧から $-\alpha^2r/\mu_0$、張力から同じ成分が出る。和は $-2\alpha^2r\hat{\vb r}/\mu_0$ である。
:::
::::

::::exercise{level=advanced}
$c_s=v_A=v$、$\theta=\pi/3$ の磁気音波の位相速度の二乗を求めよ。

:::solution
四次式を $k^4$ で割り、$u=(\omega/k)^2$ を使うと $u^2-2v^2 u+v^4/4=0$ である。従って $u=v^2(1\pm\sqrt3/2)$ となる。
:::
::::

::::exercise{level=advanced}
一様平衡で $\vb\xi=\hat{\vb y}\xi_0\cos kz$ の $\delta W/V$ を求めよ。

:::solution
$\nabla\cdot\vb\xi=0$、$\vb Q=-kB_0\xi_0\hat{\vb y}\sin kz$ である。整数波長を含む周期領域で $\sin^2$ の平均は $1/2$ なので $k^2 B_0^2\xi_0^2/(4\mu_0)$ である。
:::
::::

## 確認テスト

::::quiz
理想 MHD で物質曲面を貫く磁束の時間変化は何か。

- [x] ゼロ
- [ ] 面積が変化する場合は必ず非ゼロ
- [ ] 磁気圧が正なら非ゼロ

:::explanation
移動境界の電場 $\vb E+\vb U\times\vb B$ はゼロなので式 :eqref[eq-plasma-flux] の右辺は消える。
:::
::::

::::quiz{answer="6" tol="1e-6"}
$k_\parallel=2\,\mathrm{m^{-1}}$、$v_A=3\,\mathrm{m\,s^{-1}}$ の正の Alfvén 角周波数を答えよ。

:::explanation
式 :eqref[eq-plasma-alfven] から $\omega=k_\parallel v_A=6\,\mathrm{s^{-1}}$ である。
:::
::::

::::quiz
磁場の曲率による張力の向きはどれか。

- [x] 曲率中心へ向く。
- [ ] 常に磁場と平行である。
- [ ] 曲率中心から離れる。

:::explanation
垂直力は $B^2\vb\kappa/\mu_0$ であり、曲率ベクトルは曲率中心へ向く。
:::
::::

::::quiz{answer="5" tol="1e-6"}
$U=5\,\mathrm{m\,s^{-1}}, L=2\,\mathrm m,\eta_m=2\,\mathrm{m^2\,s^{-1}}$ の $R_m$ を答えよ。

:::explanation
定義から $UL/\eta_m=10/2=5$ である。
:::
::::
