---
title: 連続体と場の解析力学
description: 粒子列の極限から場の作用を作り、境界条件、局所保存則、正準運動方程式を導く。
notation:
  - symbol: '\mathcal L,\mathcal H'
    meaning: 空間体積あたりの Lagrangian と Hamiltonian
  - symbol: 'x^0=t'
    meaning: 秒で測る時間座標
  - symbol: '\Pi_a^\mu'
    meaning: 場の時空微分に共役な密度
  - symbol: '\pi_a'
    meaning: 場の時間微分に共役な運動量密度
---

## この章の目標

- 粒子列から弦の Lagrangian 密度と波動方程式を得る。
- 場の変分で生じる境界項から、固定端・自由端の条件を区別する。
- 場と時空の変換から Noether 流とエネルギー運動量テンソルを導く。
- 汎関数微分で場の Hamilton 方程式と等時刻の括弧を構成する。

## 弦の振動の Lagrangian

弦の中央を弾いてから離しても、振動はその場所にとどまらず端へ伝わる。
そこで弦を短い粒子列として記述し、個々の粒子のエネルギーを「単位長さあたり」に取り直す。
連続極限をとった後も、隣り合う区間の間でどのようにエネルギーが受け渡されるかを追える。

### 離散模型で次元を合わせる

長さ $l$ の一様な弦を間隔 $a=l/N$ に区切る。
平衡時の単位長さあたりの質量を $\rho$（$\mathrm{kg/m}$）、一定の張力を $\mathcal T$（$\mathrm N$）、横変位を $u_j(t)$ とする。
傾き $|u_{j+1}-u_j|/a\ll1$ で一区間の伸びは

$$
\sqrt{a^2+(u_{j+1}-u_j)^2}-a
=\frac{(u_{j+1}-u_j)^2}{2a}+O\!\left(a\left|\frac{u_{j+1}-u_j}{a}\right|^4\right).
$$

張力が既にかかっているとき、伸びに張力を掛けたものが主要な弾性エネルギーとなる。
縦運動、曲げ剛性、摩擦を無視し、固定端 $u_0=u_N=0$ とすれば

$$
L_a=\sum_{j=1}^{N-1}\frac{\rho a}2\dot u_j^2-
\sum_{j=0}^{N-1}\frac{\mathcal T}{2a}(u_{j+1}-u_j)^2 .
$$

内点での Lagrange 方程式は $\rho a\ddot u_j=\mathcal T(u_{j+1}-2u_j+u_{j-1})/a$。
空間変化の長さが $a$ より十分大きければ$u$ が空間について $C^4$ 級なら差分の分子は $a^2u_{xx}+O(a^4)$ となる。
一方、作用は Riemann 和から積分に移り

$$
L=\int_0^l\mathcal L\,\dd x,\quad
\mathcal L=\frac{\rho}{2}u_t^2-\frac{\mathcal T}{2}u_x^2,\qquad
\rho u_{tt}=\mathcal T u_{xx}
\label{eq-am-string}
$$

となる。
$\mathcal L$ を :term[Lagrangian 密度]{en="Lagrangian density" yomi="らぐらんじあんみつど"} という。
弦での単位は $\mathrm{J/m}$、三次元の体積密度なら $\mathrm{J/m^3}$ である。

### 波が伝わる速さと近似の範囲

進行波 $u=f(x-ct)$ を入れると $c=\sqrt{\mathcal T/\rho}$ である。
固定端では $\sin(n\pi x/l)$ を空間依存に取ると $\omega_n=n\pi c/l$、$n=1,2,\ldots$。
例えば $\rho=0.010\,\mathrm{kg/m}$、$\mathcal T=100\,\mathrm N$、$l=0.50\,\mathrm m$ なら $c=100\,\mathrm{m/s}$、
基本周波数は $c/(2l)=100\,\mathrm{Hz}$ である。
連続体にも基準振動があり、:ref[thm-am-normal-modes] の有限個の固有値が、弦では可算個の固有値に替わる。

離散系のまま $u_j=A\cos(kja-\omega t)$ を入れると

$$
\omega^2=\frac{4\mathcal T}{\rho a^2}\sin^2(ka/2).
$$

$|ka|\ll1$ のときは $\omega\simeq c|k|$ だが、短波長では位相速度が変わる。
例えば $ka=\pi$ では $\omega/(c|k|)=2/\pi$。
傾きが小さいという振幅の条件と、連続極限を使えるという波長の条件は別々である。

### 局所的な力と全エネルギーを別々に確かめる

内点 $j$ の運動方程式には、ばねに相当する区間が左と右から一つずつ寄与する。
連続極限で一次の張力そのものが加わるのでなく、その空間差が質量を加速する。
傾きが空間的に一定なら弦は傾いていても二つの張力の横成分が相殺し、加速度は零である。
山の頂上では $u_{xx}<0$ であり、張力は頂上を下向きに加速する。

:::example[短い区間を通過するパルス]
無限長の一様弦に $u(x,t)=A\exp[-(x-ct)^2/(2d^2)]$ を置く。
$u_t=-cu_x$ だから波動方程式を満たす。
$x=ct$ では変位は最大だが速度も傾きも零で、そこだけを見ればエネルギー密度は零である。
一方、最大傾斜の位置は $|x-ct|=d$ で、$|u_x|_{\max}=(A/d)e^{-1/2}$。
$A/d=0.02$ なら最大傾斜は約 $0.012$ で、長さの展開近似は首尾一貫している。
パルスを「山の高さだけでエネルギーが決まる」と考えてはならない。
:::

### 端に達した波が返る向き

$x\le0$ の弦で端が $x=0$ にあるとする。
右へ進む入射波を $f(t-x/c)$、左へ戻る波を $g(t+x/c)$ とする。
端 $x=0$ が固定なら $f(t)+g(t)=0$ なので $g=-f$。
変位の符号が反転する。
自由端なら $u_x(0,t)=[-f'(t)+g'(t)]/c=0$ なので $g=f+\mathrm{const}$ となり、静的な変位の原点を合わせれば符号を保って反射する。
この違いは次節の作用の端点項だけからも決められる。

### 初期変位・初速度から波を組み立てる

無限に長い一様弦で波動方程式を初期値問題として解く。特性座標 $\xi=x-ct,\eta=x+ct$ に移ると $\partial_t^2-c^2\partial_x^2=-4c^2\partial_\xi\partial_\eta$。従って滑らかな解は $u=F(\xi)+G(\eta)$ である。初期変位 $f(x)$ と初速度 $g(x)$ を使うと $F+G=f,\ c(-F'+G')=g$。第一式を微分して第二式との和差を取れば $F'=(f'-g/c)/2,\ G'=(f'+g/c)/2$ となる。積分定数は $F+G=f$ により相殺し、
$$
u(x,t)=\tfrac12[f(x-ct)+f(x+ct)]+\frac1{2c}\int_{x-ct}^{x+ct}g(y)\,\dd y
\label{eq-am-string-cauchy}
$$
を得る。波の到着を決めるのは初期変位では二端の値、初速度では区間全体である。初め変形のない弦の短い区間だけをたたくと、パルスの中間に一定変位が残ることがある。その部分に横速度も傾きもなければエネルギー密度は零である。

:::example[定在波を静止状態から放す]
弦上の位置は m、時刻は s とし、三角関数の引数はこれらの数値で記す。線密度 $\rho=1\,\mathrm{kg/m}$、波速 $c=2\,\mathrm{m/s}$、初期変位 $f=\frac{\cos{\left(x \right)}}{100}\,\mathrm m$、初速度 $g=0\,\mathrm{m/s}$ とする。二つの特性の寄与を初期値から決めれば

$$
u(x,t)=- \frac{\left(2 \sin^{2}{\left(t \right)} - 1\right) \cos{\left(x \right)}}{100},\quad u_t=- \frac{\sin{\left(t \right)} \cos{\left(t \right)} \cos{\left(x \right)}}{25}.
$$

$(x,t)=(\frac{\pi}{2},\frac{\pi}{8})$ での変位、速度、傾きは順に $0,0,- \frac{\sqrt{2}}{200}$。この位置でのエネルギー密度は $(\rho u_t^2+\mathcal T u_x^2)/2=\frac{1}{10000}\,\mathrm{J/m}$、流束は $-\mathcal T u_tu_x=0\,\mathrm W$。速度を与えずに離すときは、二つの等振幅の進行波が必要である。観測点が節なら変位が常に零でも、傾きによる弾性エネルギーをもちうる。
:::

:::example[変位ゼロの弦に速度の分布を与える]
弦上の位置は m、時刻は s とし、三角関数の引数はこれらの数値で記す。線密度 $\rho=1\,\mathrm{kg/m}$、波速 $c=3\,\mathrm{m/s}$、初期変位 $f=0\,\mathrm m$、初速度 $g=\frac{3 \cos{\left(2 x \right)}}{100}\,\mathrm{m/s}$ とする。二つの特性の寄与を初期値から決めれば

$$
u(x,t)=\frac{\sin{\left(6 t - 2 x \right)} + \sin{\left(6 t + 2 x \right)}}{400},\quad u_t=\frac{3 \left(\cos{\left(6 t - 2 x \right)} + \cos{\left(6 t + 2 x \right)}\right)}{200}.
$$

$(x,t)=(0,\frac{\pi}{12})$ での変位、速度、傾きは順に $\frac{1}{200},0,0$。この位置でのエネルギー密度は $(\rho u_t^2+\mathcal T u_x^2)/2=0\,\mathrm{J/m}$、流束は $-\mathcal T u_tu_x=0\,\mathrm W$。初速度の積分が正の範囲で変位が育つ。振動の四分の一周期で運動から弾性変形へエネルギーが移る。
:::

:::example[一方向だけへ進む波を準備する]
弦上の位置は m、時刻は s とし、三角関数の引数はこれらの数値で記す。線密度 $\rho=1\,\mathrm{kg/m}$、波速 $c=2\,\mathrm{m/s}$、初期変位 $f=\frac{\sin{\left(x \right)}}{100}\,\mathrm m$、初速度 $g=- \frac{\cos{\left(x \right)}}{50}\,\mathrm{m/s}$ とする。二つの特性の寄与を初期値から決めれば

$$
u(x,t)=- \frac{\sin{\left(2 t - x \right)}}{100},\quad u_t=- \frac{\cos{\left(2 t - x \right)}}{50}.
$$

$(x,t)=(\frac{\pi}{2},\frac{\pi}{4})$ での変位、速度、傾きは順に $0,- \frac{1}{50},\frac{1}{100}$。この位置でのエネルギー密度は $(\rho u_t^2+\mathcal T u_x^2)/2=\frac{1}{2500}\,\mathrm{J/m}$、流束は $-\mathcal T u_tu_x=\frac{1}{1250}\,\mathrm W$。初速度を $-cf'$ に合わせると左へ動く成分が消える。この準備がない限り、与えた形をそのまま一方向へ平行移動した式は初期条件を満たさない。
:::

:::example[別々の波長を変位と速度に入れる]
弦上の位置は m、時刻は s とし、三角関数の引数はこれらの数値で記す。線密度 $\rho=1\,\mathrm{kg/m}$、波速 $c=1\,\mathrm{m/s}$、初期変位 $f=\frac{\sin{\left(x \right)}}{100}\,\mathrm m$、初速度 $g=\frac{\cos{\left(2 x \right)}}{50}\,\mathrm{m/s}$ とする。二つの特性の寄与を初期値から決めれば

$$
u(x,t)=\frac{- \sin{\left(t - x \right)} + \sin{\left(t + x \right)} + \sin{\left(2 t - 2 x \right)} + \sin{\left(2 t + 2 x \right)}}{200},\quad u_t=\frac{- \cos{\left(t - x \right)} + \cos{\left(t + x \right)} + 2 \cos{\left(2 t - 2 x \right)} + 2 \cos{\left(2 t + 2 x \right)}}{200}.
$$

$(x,t)=(0,\frac{\pi}{4})$ での変位、速度、傾きは順に $\frac{1}{100},0,\frac{\sqrt{2}}{200}$。この位置でのエネルギー密度は $(\rho u_t^2+\mathcal T u_x^2)/2=\frac{1}{40000}\,\mathrm{J/m}$、流束は $-\mathcal T u_tu_x=0\,\mathrm W$。変位と速度の二つの Fourier 成分は独立に進み、エネルギー密度には干渉項が現れる。一つの点の変位だけから全波形のエネルギーを逆算できない。
:::

:eqref[eq-am-string-cauchy] の積分区間を時空図に表すと、観測点から過去に向かう二本の特性の足にはさまれた部分になる。実際の弦に端があり、三角形が端に達する時刻を越えた場合は、初期線だけでは値を決められない。固定端なら端を越える初期値を奇関数、自由端なら偶関数として延長すれば、反射条件を満たしながら同じ図を利用できる。

```tikz caption="初期値問題の依存領域。二本の特性で切り取った初期線上の区間だけが観測点に影響する。"
\draw[->] (-2.5,0)--(2.7,0) node[right] {$x$};
\draw[->] (-2.35,-.15)--(-2.35,2.6) node[above] {$t$};
\fill[gray!15] (-1.6,0)--(1.9,0)--(.15,2.05)--cycle;
\draw[thick,->] (-1.6,0)--(.15,2.05);
\draw[thick,->] (1.9,0)--(.15,2.05);
\fill (.15,2.05) circle (1.6pt) node[above] {$(x,t)$};
\node[below] at (-1.6,0) {$x-ct$};
\node[below] at (1.9,0) {$x+ct$};
```

## 場の Euler–Lagrange 方程式

```tikz caption="弦の区間が受ける張力の向き。右端での横成分から左端での横成分を引くと、区間の運動量の時間変化になる。"
\draw[->] (-.25,0) -- (5.1,0) node[right] {$x$};
\draw[thick] (0,1.55) .. controls (1,.95) and (1.6,.75) .. (2,.8)
 .. controls (3,.9) and (3.6,1.2) .. (4.6,1.85);
\draw[dashed] (1,0) -- (1,1.08);
\draw[dashed] (3.65,0) -- (3.65,1.32);
\draw[->,thick] (1,1.08) -- (.25,1.38);
\draw[->,thick] (3.65,1.32) -- (4.4,1.74);
\node[below] at (1,0) {$a$};
\node[below] at (3.65,0) {$b$};
\node[above] at (.53,1.45) {$\mathcal T$};
\node[above] at (4.1,1.7) {$\mathcal T$};
```

### 場の値を一つの関数だけ動かす

空間領域 $D\subset\R^d$ の各点で実数 $\phi^a(\vb x,t)$ を指定する。
$a=1,\ldots,s$ は場の成分で、$\mu=0,\ldots,d$ は独立変数の番号とする。
$x^0=t$、$\partial_0=\partial_t$ なので時間の単位は秒である。
Lorentz 指標の上げ下げは用いない。
積分要素を $\dd X=\dd t\,\dd^dx$ とし、一次微分までの局所作用を

$$
S[\phi]=\int_{t_1}^{t_2}\!\int_D\mathcal L(\phi^a,\partial_\mu\phi^a,X)\,\dd^dx\,\dd t,\qquad
\Pi_a^\mu=\frac{\partial\mathcal L}{\partial(\partial_\mu\phi^a)}
$$

とおく。
$\mu$ が添えられた引数は全成分の組を表し、縮約する場合は和を明示する。

$\phi^a\mapsto\phi^a+\epsilon\eta^a$ とし、時間端では $\eta^a=0$ とする。
第一変分で $\partial_\mu\eta^a$ を一度部分積分すると

$$
\begin{aligned}
\delta S={}&\int\sum_{a=1}^s
\left(\mathcal L_{\phi^a}-\sum_{\mu=0}^d\partial_\mu\Pi_a^\mu\right)\eta^a\,\dd X\\
&+\int_{t_1}^{t_2}\!\int_{\partial D}
\sum_{a=1}^s\sum_{i=1}^d n_i\Pi_a^i\eta^a\,\dd A\,\dd t .
\end{aligned}
$$

内部に台を持つ任意の滑らかな $\eta^a$ を使うと、括弧内が連続で一点で正なら、その近くに非負の試験関数を置いて積分を正にできる。
負の場合も符号を逆にすればよい。
これは :ref[lem-am-fundamental-variation] の多変数版である。
従って局所方程式は

$$
\mathcal E_a\equiv\mathcal L_{\phi^a}-\sum_{\mu=0}^d\partial_\mu\Pi_a^\mu=0 .
\label{eq-am-field-el}
$$

### 境界を選ぶと何が変わるか

場を境界で指定すると $\eta^a=0$。
場を自由にする境界で外部作用がなければ $\sum_i n_i\Pi_a^i=0$ が自然境界条件である。
弦の右端なら $\Pi^x=-\mathcal T u_x$ より自由端で $u_x=0$。
周期条件では両端の場と流束が一致し、境界項が差として消える。

:::example[端をばねに結ぶ]
右端にばね定数 $\kappa$ のばねを取り付け、$L$ に $-\kappa u(l)^2/2$ を加える。
右端での変分は $[-\mathcal T u_x(l)-\kappa u(l)]\eta(l)$ だから $\mathcal T u_x(l)+\kappa u(l)=0$。
例えば $\mathcal T=20\,\mathrm N$、$\kappa=100\,\mathrm{N/m}$、$u(l)=0.01\,\mathrm m$ なら端の傾きは $-0.05$ である。
端のばね力と張力の横成分がつり合う。
:::

:::example[局所的な復元力]
床ばねをもつ弦では $\mathcal L=\rho u_t^2/2-\mathcal T u_x^2/2-\kappa u^2/2$ とする。
この $\kappa$ の単位は $\mathrm{N/m^2}$ であり、端ばねの定数と異なる。
局所方程式は $\rho u_{tt}=\mathcal T u_{xx}-\kappa u$、平面波では $\omega^2=(\mathcal T k^2+\kappa)/\rho$ である。
$k=0$ でも復元力が残る。
:::

### 材料の界面で必要な条件

二本の弦を $x=0$ で結ぶ。張力が左右で異なる場合、結び目を横方向には自由に動かせる案内によって、長手方向の力だけを支えるとする。
結合が外れないので変位は連続であり、結び目の許容変分 $\eta(0)$ は左右で同じである。結び目の集中質量と横方向の外力は零とする。
両側を別々に部分積分すると左の右端と右の左端は逆向きになるため

$$
[-\mathcal T_1u_x(0^-)+\mathcal T_2u_x(0^+)]\eta(0)=0 .
$$

従って $\mathcal T_1u_x(0^-)=\mathcal T_2u_x(0^+)$。
傾きそのものが等しいのは張力も等しいときである。

:::example[二本の弦の反射率と透過率]
正弦波の複素表示で、左は $u_1=(e^{ik_1x}+r e^{-ik_1x})e^{-i\omega t}$、右は $u_2=t_a e^{ik_2x-i\omega t}$ とする。
$c_j=\sqrt{\mathcal T_j/\rho_j}$、$k_j=\omega/c_j$ で、
$Z_j=\mathcal T_j/c_j=\sqrt{\mathcal T_j\rho_j}$ を波動インピーダンスという。
変位と張力の連続条件は

$$
1+r=t_a,\qquad Z_1(1-r)=Z_2t_a .
$$

よって $r=(Z_1-Z_2)/(Z_1+Z_2)$、$t_a=2Z_1/(Z_1+Z_2)$。
$Z_2=3Z_1$ なら $r=-1/2,\ t_a=1/2$。
時間平均のエネルギー流は $Z\omega^2|A|^2/2$ に比例するので、入射に対する反射は $1/4$、透過は $3/4$ となる。
透過エネルギーの割合を単に $|t_a|^2$ とすると、異なる二本の弦のインピーダンスを取り落とす。
:::

### 全発散を加えた密度

任意の滑らかな $C^\mu(\phi,X)$ に対し $\mathcal L'=\mathcal L+\sum_\mu\partial_\mu C^\mu$ とする。
作用差は時空領域の境界積分である。
境界で場を固定した変分ではこの差の変分は零なので、内部の Euler–Lagrange 方程式は同じである。
一方、外力を与える自然境界条件まで含める問題では、境界作用の変化を同時に数える必要がある。

例として弦で $\mathcal L'=\mathcal L+\partial_x(\kappa u^2/2)$ とすると、追加項は $\kappa uu_x$。
その Euler–Lagrange 項は $\kappa u_x-\partial_x(\kappa u)=0$（$\kappa$ 一定）であるが、右端での変分に $\kappa u(l)\eta(l)$ が加わる。
従って同じ自由端を記述するには、対応する端点項も差し引かなければならない。

### 多成分の境界条件を行列で書く

二種の場の勾配がエネルギーで結合する場合、自然境界条件が各成分の法線微分零になるとは限らない。
一方向の二成分弦で、列 $u=(u_1,u_2)^\top$、一定の対称正定値行列 $M,K$ に対し

$$
\mathcal L=\tfrac12 u_t^\top Mu_t-\tfrac12u_x^\top Ku_x
$$

とする。
$\Pi^0=Mu_t,\ \Pi^x=-Ku_x$ なので局所方程式は $Mu_{tt}=Ku_{xx}$。
右端で二成分とも独立に自由なら $Ku_x=0$ であり、$K$ が可逆なので $u_x=0$ である。
一成分だけ、例えば $u_1$ を固定し $u_2$ が自由なら、$\delta u_1=0$ と $(Ku_x)_2=0$ だけを課す。
二成分の境界力を別々に零と置くと過剰な条件になる。

例えば $K=\mathcal T\begin{pmatrix}2&1\\1&2\end{pmatrix}$ とすれば第二条件は $u_{1,x}+2u_{2,x}=0$ となる。
$u_1$ を固定した支持具は力を担えるので、そこで $(Ku_x)_1$ が零である必要はない。
拘束で許された仮想変位に対する仕事だけを零とする :ref[def-am-ideal-constraint] と同じ考え方である。

### 等方弾性体の変位場と境界に働く力

固体の各点の三次元変位を $u_i(\vb x,t)$ とし、基準配置の座標 $x^i$ を m で測る。小さい変位勾配の一次まででは、微小線分の長さの二乗の変化は $2\sum_{i,j}e_{ij}\dd x^i\dd x^j$、$e_{ij}=(\partial_i u_j+\partial_j u_i)/2$ となる。変位勾配の反対称部分はその点での微小回転を表し、線分の長さの一次変化には寄与しない。実対称な $e$ を直交対角化し、等方性を課すとその二次エネルギーは三つの固有値の対称な二次式になる。従って $(\tr e)^2$ と $\tr(e^2)$ の線形結合しかない。体積あたりの位置エネルギーを
$$
W=\frac\lambda2\left(\sum_i e_{ii}\right)^2+\mu\sum_{i,j}e_{ij}^2
$$
とする。$\lambda,\mu$ は Lamé 定数で Pa の単位をもつ。ここで $\lambda$ は拘束の乗数とは別の材料定数である。$e$ を球状部分と跡のない部分に分けると、エネルギーの正定値条件は $\mu>0,\ K=\lambda+2\mu/3>0$ である。$K$ は体積弾性率であり、$\lambda$ 単独に正を要求する必要はない。

$\mathcal L=\rho\sum_i\dot u_i^2/2-W+\sum_i f_i u_i$ として $\delta W=\sum_{i,j}\sigma_{ij}\delta e_{ij}$ と書くと
$$
\sigma_{ij}=\lambda\left(\sum_k e_{kk}\right)\delta_{ij}+2\mu e_{ij},\qquad
\rho\ddot u_i=\sum_j\partial_j\sigma_{ij}+f_i .
\label{eq-am-elastic-stress}
$$
実際、$\sigma$ の対称性によって $\delta W=\sum_{i,j}\sigma_{ij}\partial_j\delta u_i$。これを空間で部分積分すると内部に上の運動方程式、境界に $-\sum_{i,j}\sigma_{ij}n_j\delta u_i$ が残る。表面力 $\vb t$ を指定するなら境界に $\int\vb t\cdot\vb u\,\dd S$ を加え、$\sum_j\sigma_{ij}n_j=t_i$。面の外向き法線を変えると力の向きも反転する。応力の非対角成分はせん断を測り、対称性は微小体積の角運動量の収支とも整合する。

:::example[体積だけを増やす膨張]
材料定数を $\lambda=2.0\,\mathrm{GPa},\mu=1.0\,\mathrm{GPa}$ とする。小さいひずみと切断面の外向き法線を

$$
e=\left(\begin{matrix}\frac{1}{1000} & 0 & 0\\0 & \frac{1}{1000} & 0\\0 & 0 & \frac{1}{1000}\end{matrix}\right),\qquad \vb n=\left(\begin{matrix}1\\0\\0\end{matrix}\right).
$$

とする。このとき体積の相対変化は一次まで $\tr e=\frac{3}{1000}$ である。応力の球状項を先に作って $2\mu e$ を足せば

$$
\sigma=\left(\begin{matrix}8 & 0 & 0\\0 & 8 & 0\\0 & 0 & 8\end{matrix}\right)\,\mathrm{MPa},\qquad \vb t=\sigma\vb n=\left(\begin{matrix}8\\0\\0\end{matrix}\right)\,\mathrm{MPa}.
$$

法線成分は $\vb n\cdot\vb t=8\,\mathrm{MPa}$。単位体積の蓄積エネルギーは $W=\sum_{i,j}\sigma_{ij}e_{ij}/2=12000\,\mathrm{J/m^3}$ であり、段階的にひずみを零から増やしたときの仕事と一致する。三方向の伸びが等しいので表面力は各面の法線に平行である。圧力を圧縮で正と定義する場合はこの引張応力と符号が逆になる。
:::

:::example[体積を変えない単純せん断]
材料定数を $\lambda=2.0\,\mathrm{GPa},\mu=1.0\,\mathrm{GPa}$ とする。小さいひずみと切断面の外向き法線を

$$
e=\left(\begin{matrix}0 & \frac{1}{1000} & 0\\\frac{1}{1000} & 0 & 0\\0 & 0 & 0\end{matrix}\right),\qquad \vb n=\left(\begin{matrix}0\\1\\0\end{matrix}\right).
$$

とする。このとき体積の相対変化は一次まで $\tr e=0$ である。応力の球状項を先に作って $2\mu e$ を足せば

$$
\sigma=\left(\begin{matrix}0 & 2 & 0\\2 & 0 & 0\\0 & 0 & 0\end{matrix}\right)\,\mathrm{MPa},\qquad \vb t=\sigma\vb n=\left(\begin{matrix}2\\0\\0\end{matrix}\right)\,\mathrm{MPa}.
$$

法線成分は $\vb n\cdot\vb t=0\,\mathrm{MPa}$。単位体積の蓄積エネルギーは $W=\sum_{i,j}\sigma_{ij}e_{ij}/2=2000\,\mathrm{J/m^3}$ であり、段階的にひずみを零から増やしたときの仕事と一致する。工学せん断ひずみは $2e_{xy}$ である。非対角成分を二回数える二重和とこの二倍を混同するとエネルギーが合わない。
:::

:::example[横方向を固定した一軸伸長]
材料定数を $\lambda=2.0\,\mathrm{GPa},\mu=1.0\,\mathrm{GPa}$ とする。小さいひずみと切断面の外向き法線を

$$
e=\left(\begin{matrix}\frac{1}{1000} & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{matrix}\right),\qquad \vb n=\left(\begin{matrix}0\\1\\0\end{matrix}\right).
$$

とする。このとき体積の相対変化は一次まで $\tr e=\frac{1}{1000}$ である。応力の球状項を先に作って $2\mu e$ を足せば

$$
\sigma=\left(\begin{matrix}4 & 0 & 0\\0 & 2 & 0\\0 & 0 & 2\end{matrix}\right)\,\mathrm{MPa},\qquad \vb t=\sigma\vb n=\left(\begin{matrix}0\\2\\0\end{matrix}\right)\,\mathrm{MPa}.
$$

法線成分は $\vb n\cdot\vb t=2\,\mathrm{MPa}$。単位体積の蓄積エネルギーは $W=\sum_{i,j}\sigma_{ij}e_{ij}/2=2000\,\mathrm{J/m^3}$ であり、段階的にひずみを零から増やしたときの仕事と一致する。横変位を固定するため、横の面にも垂直な反力が必要である。横の面を自由にした引張試験とは違う条件である。
:::

:::example[せん断と膨張を受ける斜面]
材料定数を $\lambda=2.0\,\mathrm{GPa},\mu=1.0\,\mathrm{GPa}$ とする。小さいひずみと切断面の外向き法線を

$$
e=\left(\begin{matrix}\frac{1}{500} & \frac{1}{1000} & 0\\\frac{1}{1000} & - \frac{1}{1000} & 0\\0 & 0 & 0\end{matrix}\right),\qquad \vb n=\left(\begin{matrix}\frac{3}{5}\\\frac{4}{5}\\0\end{matrix}\right).
$$

とする。このとき体積の相対変化は一次まで $\tr e=\frac{1}{1000}$ である。応力の球状項を先に作って $2\mu e$ を足せば

$$
\sigma=\left(\begin{matrix}6 & 2 & 0\\2 & 0 & 0\\0 & 0 & 2\end{matrix}\right)\,\mathrm{MPa},\qquad \vb t=\sigma\vb n=\left(\begin{matrix}\frac{26}{5}\\\frac{6}{5}\\0\end{matrix}\right)\,\mathrm{MPa}.
$$

法線成分は $\vb n\cdot\vb t=\frac{102}{25}\,\mathrm{MPa}$。単位体積の蓄積エネルギーは $W=\sum_{i,j}\sigma_{ij}e_{ij}/2=8000\,\mathrm{J/m^3}$ であり、段階的にひずみを零から増やしたときの仕事と一致する。面の向きを指定すると同じ応力テンソルから法線成分と接線成分が同時に決まる。
:::

### 変位の偏光と二つの音速

材料定数が空間的に一定で体積力がないとき、応力の発散を展開すると $\rho\ddot{\vb u}=\mu\nabla^2\vb u+(\lambda+\mu)\nabla(\nabla\cdot\vb u)$。平面波 $\vb u=\vb A e^{i(\vb k\cdot\vb x-\omega t)}$ の実部を取れば
$$
\rho\omega^2\vb A=\mu k^2\vb A+(\lambda+\mu)\vb k(\vb k\cdot\vb A).
$$
$\vb A$ を $\vb k$ に平行と垂直に分解すると、縦波が一つで $c_L=\sqrt{(\lambda+2\mu)/\rho}$、横波が二つで $c_T=\sqrt{\mu/\rho}$。材料の回転対称性により二つの横偏光の速度は等しい。例えば $\lambda=2\,\mathrm{GPa},\mu=1\,\mathrm{GPa},\rho=1000\,\mathrm{kg/m^3}$ とすると $c_L=2000,\ c_T=1000\,\mathrm{m/s}$。同じ時刻に励起し、界面や反射のない厚さ $0.10\,\mathrm m$ の媒質を進むなら最初の到着時間は $50$ と $100\,\mathrm{\mu s}$ である。この議論は体積中を進む波に限られ、自由表面に沿う波には境界条件を加えて解く必要がある。

### 端に慣性のある弦の固有値と規格化

端質量があると通常の弦の内積でモードを直交させることはできない。左固定、右端に質量 $M$ とばね $\kappa$ があるとし、時間に依存しない線密度 $\rho(x)>0$、張力 $\mathcal T(x)>0$ を取る。$u(x,t)=f(x)\cos\omega t$ の内部で $-(\mathcal Tf')'=\omega^2\rho f$。空間微分を含む式は張力が位置に依存する場合にも使える。端の質量にも Lagrange 方程式を適用すると
$$
f(0)=0,\qquad \mathcal T(l) f'(l)+\kappa f(l)=M\omega^2 f(l).
$$
従って質量効果を端ばねと同じ正の剛性として足してはならない。周波数が高いほど慣性力の符号と大きさが効く。異なる固有値のモード $f_n,f_m$ の方程式に相手を掛け、引いて区間全体で積分する。内部の $\int\mathcal T f_n'f_m'$ は対称なので消え、右端での二項に上の条件を用いると
$$
(\omega_n^2-\omega_m^2)
\left[\int_0^l \rho f_n f_m\,\dd x+M f_n(l)f_m(l)\right]=0 .
$$
これは端でも横方向の運動エネルギーを蓄えることの直接の結果である。$M=0$ と置けば分布した質量だけの内積になる。今度は単一モードで $f$ 自身を掛けて部分積分すれば
$$
\omega^2=\frac{\int_0^l\mathcal T(f')^2\,\dd x+\kappa f(l)^2}
{\int_0^l\rho f^2\,\dd x+M f(l)^2}.
\label{eq-am-string-end-mass}
$$
端のばねと端の質量はそれぞれ分子、分母に入る。分母の内積でモードを規格化して変位と初速度を射影すると、端の初速度からもモードが励起される。

一様弦で $\kappa=0$ なら $f=\sin(kx)$、$z=kl$ として $\cot z=\eta z,\ \eta=M/(\rho l)$。左辺は各区間 $(n\pi,n\pi+\pi/2)$ で正の無限大から零まで単調減少し、右辺は増加するので $\eta>0$ では各区間に一つの根がある。その右側の半区間では符号が逆で根はない。形式的に捨てた $z=0$ は $f=Cx$ で境界条件から $C=0$ となり、固有値には含まれない。$\eta\gg1$ の基本根は $z\ll1$ であり、$1/z-z/3+O(z^3)=\eta z$ より $z^2\simeq1/(\eta+1/3)$。従って基本角振動数の二乗は $\mathcal T/[l(M+\rho l/3)]$ に近い。弦がほぼ直線状に動けば、速度は端の値の $x/l$ 倍、運動エネルギーを端速度で測った有効質量は $\int_0^l\rho(x/l)^2\dd x=\rho l/3$ である。この一致により位相と慣性の双方から近似の範囲を判断できる。

:::example[軽い端質量]
$l=1\,\mathrm m,\rho=0.10\,\mathrm{kg/m},\mathcal T=10\,\mathrm N$、右端の質量を $M=\frac{1}{100}\,\mathrm{kg}$ とする。$c=10\,\mathrm{m/s}$、$\eta=\frac{1}{10}$ である。端の式を $\cos z-\eta z\sin z=0$ とし、上で示した区間で二分すれば基本根と次の根は $1.428870,\ 4.305801$。それぞれに $c/l$ を掛けた $14.28870,\ 43.05801\,\mathrm{s^{-1}}$ が角振動数である。直線状の試行関数を商に入れて求めた基本角振動数は $15.19109\,\mathrm{s^{-1}}$。厳密な基本根に対する相対差は約 $6.3154\%$ であり、質量をもつ端の動きと弦の変形を一つの直線で表したときの精度を見積もれる。
:::

:::example[弦と等しい端質量]
$l=1\,\mathrm m,\rho=0.10\,\mathrm{kg/m},\mathcal T=10\,\mathrm N$、右端の質量を $M=\frac{1}{10}\,\mathrm{kg}$ とする。$c=10\,\mathrm{m/s}$、$\eta=1$ である。端の式を $\cos z-\eta z\sin z=0$ とし、上で示した区間で二分すれば基本根と次の根は $0.860334,\ 3.425618$。それぞれに $c/l$ を掛けた $8.60334,\ 34.25618\,\mathrm{s^{-1}}$ が角振動数である。直線状の試行関数を商に入れて求めた基本角振動数は $8.66025\,\mathrm{s^{-1}}$。厳密な基本根に対する相対差は約 $0.6616\%$ であり、質量をもつ端の動きと弦の変形を一つの直線で表したときの精度を見積もれる。
:::

:::example[重い端質量]
$l=1\,\mathrm m,\rho=0.10\,\mathrm{kg/m},\mathcal T=10\,\mathrm N$、右端の質量を $M=1\,\mathrm{kg}$ とする。$c=10\,\mathrm{m/s}$、$\eta=10$ である。端の式を $\cos z-\eta z\sin z=0$ とし、上で示した区間で二分すれば基本根と次の根は $0.311053,\ 3.173097$。それぞれに $c/l$ を掛けた $3.11053,\ 31.73097\,\mathrm{s^{-1}}$ が角振動数である。直線状の試行関数を商に入れて求めた基本角振動数は $3.11086\,\mathrm{s^{-1}}$。厳密な基本根に対する相対差は約 $0.0105\%$ であり、質量をもつ端の動きと弦の変形を一つの直線で表したときの精度を見積もれる。
:::

### 曲げを受ける細いはりでは何を端で固定するか

弦は傾きにエネルギーを蓄える。一方、中立軸の曲率に比例した曲げモーメントをもつ細いはりでは、小さい横変位の曲率は $u_{xx}$ である。長さに比べて断面が小さく、せん断変形と断面の回転慣性を無視できる低周波を扱う。線密度 $\rho$、曲げ剛性 $B$（$\mathrm{N\,m^2}$）、単位長さの外力 $f$ に対して $\mathcal L=\rho u_t^2/2-Bu_{xx}^2/2+fu$ とする。空間で二回部分積分すれば
$$
\delta S=\int\!\!\int[-\rho u_{tt}-(Bu_{xx})_{xx}+f]\delta u\,\dd x\,\dd t
+\int[(Bu_{xx})_x\delta u-Bu_{xx}\delta u_x]_0^l\,\dd t .
$$
従って Euler–Bernoulli 方程式は $\rho u_{tt}+(Bu_{xx})_{xx}=f$。固定端では変位と傾きの両方を固定する。単純支持端なら変位は固定だが傾きは自由なので $u=0,\ Bu_{xx}=0$。自由端なら二つの変分が独立で $Bu_{xx}=0,\ (Bu_{xx})_x=0$。それぞれ曲げモーメントがない条件とせん断力がない条件である。弦の自由端の $u_x=0$ をそのまま使うことはできない。

一様な両端単純支持はりの $\sin(n\pi x/l)$ は変位と二階微分が両端で消える。代入すれば $\omega_n=\sqrt{B/\rho}(n\pi/l)^2$、第二モードの周波数は基本の四倍である。例えば $B=0.25\,\mathrm{N\,m^2},\rho=1\,\mathrm{kg/m},l=1\,\mathrm m$ なら $\omega_1=\pi^2/2\simeq4.93\,\mathrm{s^{-1}}$、基本周波数は約 $0.785\,\mathrm{Hz}$。同じ長さでも弦はモード番号の一乗で増えるから、振動数比を測るだけで曲げ支配と張力支配を見分ける手がかりになる。両方がある場合には、分散関係の右辺に張力由来の $k^2$ と曲げ由来の $k^4$ が加わる。

## 場の Noether の定理

### 座標セルの体積変化も含める

粒子の :ref[thm-am-noether] では時間区間の伸びが必要だった。
場の場合には時間と空間の両方が変換され、積分領域そのものも変わる。
微小変換を $x'^\mu=x^\mu+\epsilon\xi^\mu(X)$、対応する点で $\phi'^a(X')=\phi^a(X)+\epsilon\Delta^a(\phi,X)$ とする。
同じ $X$ での変化の係数は

$$
\eta^a=\Delta^a-\sum_{\nu=0}^d\xi^\nu\partial_\nu\phi^a .
$$

微分の変化は $\partial_\mu\Delta^a-\sum_\nu(\partial_\mu\xi^\nu)\partial_\nu\phi^a$、体積の変化率は $\sum_\mu\partial_\mu\xi^\mu$。
前者は逆 Jacobi 行列を一次展開し、後者は行列式を一次展開すれば得られる。

:::theorem[場の Noether 流]{#thm-am-field-noether}
作用密度 $\mathcal L\,\dd X$ の変換の一次変化が、任意の場について $\epsilon\sum_\mu\partial_\mu B^\mu\,\dd X$ であるとする。
すると

$$
j^\mu=\sum_{a=1}^s\Pi_a^\mu\eta^a+\mathcal L\xi^\mu-B^\mu,\qquad
\partial_t j^0+\sum_{i=1}^d\partial_i j^i=0
\label{eq-am-field-current}
$$

が解の上で成り立つ。
$j^0$ は保存量の密度、$j^i$ は流束であり、合わせて :term[Noether 流]{en="Noether current" yomi="ねーたーりゅう"} という。
:::

:::proof
固定した点での密度の変化に、動くセルが運ぶ密度の発散を加える。
座標の陽な依存もセルの項に含めれば、変換の一次係数は

$$
\sum_a\mathcal L_{\phi^a}\eta^a+
\sum_{a,\mu}\Pi_a^\mu\partial_\mu\eta^a+
\sum_\mu\partial_\mu(\mathcal L\xi^\mu).
$$

第二項を積の微分と $-\sum_{a,\mu}(\partial_\mu\Pi_a^\mu)\eta^a$ に分け、仮定の発散を引くと、

$$
0=\sum_a\mathcal E_a\eta^a+\sum_{\mu=0}^d\partial_\mu j^\mu .
$$

従って :eqref[eq-am-field-el] を満たす場では保存則を得る。
固定領域で積分すれば $\dd(\int_D j^0\dd^dx)/\dd t=-\int_{\partial D}\sum_i n_i j^i\dd A$。
全量の保存には境界を横切る正味の流束が零であることも要る。
:::

:::example[変位の原点をずらす]
床ばねのない一様弦で $u'=u+\epsilon$ とする。
$\xi=0,\Delta=1,B=0$ で密度は不変である。
流れは $j^0=\rho u_t,\ j^x=-\mathcal T u_x$。
横方向の運動量は張力の横成分を通して受け渡され、連続の式は波動方程式になる。
固定端では変位の並進が境界条件を変えるため、弦だけの全横運動量は一般に保存しない。
:::

:::example[二成分の内部回転]
二つの横変位 $u,v$ が同じ $\rho,\mathcal T$ をもち、密度が両成分の弦密度の和であるとする。
$(\Delta u,\Delta v)=(-v,u)$ は断面内の回転であり、保存密度は $\rho(uv_t-vu_t)$、流束は $-\mathcal T(uv_x-vu_x)$。
例えば $u=A\cos(\omega t-kx),\ v=A\sin(\omega t-kx)$ と置くと、各点で変位は円を描き、密度は $\rho A^2\omega$ である。
両成分を波動方程式に入れると、許される周波数は $\omega^2=\mathcal T k^2/\rho$ で決まる。
:::

### 境界項を持つ内部対称性

床ばねなしの弦で $u'=u+\epsilon t$ とすると $\Delta=t$、$\xi=0$。
$\mathcal L$ の変化は $\epsilon\rho u_t=\epsilon\partial_t(\rho u)$ である。
従って $B^0=\rho u,\ B^x=0$ と選び

$$
j^0=\rho(tu_t-u),\qquad j^x=-\mathcal Tt u_x
$$

を得る。
発散を実際にとると $t(\rho u_{tt}-\mathcal T u_{xx})$ である。
周期境界の弦では全量が保存し、空間平均の変位が等速に移動することに対応する。
有限個の自由粒子の Galilei 変換で $mx$ の境界項が現れたのと同じ仕組みである。

:::example[保存流から局所的な力の符号を確かめる]
弦の区間 $[a,b]$ の全横運動量を $P_\perp=\int_a^b\rho u_t\dd x$ とする。
流束を用いれば

$$
\dot P_\perp=\mathcal T[u_x(b)-u_x(a)] .
$$

右端での正の傾きは区間を上へ引き、左端での正の傾きは区間を下へ引く。
$\mathcal T=50\,\mathrm N$、右端の傾き $0.04$、左端 $-0.02$ なら合力は $3\,\mathrm N$。
場の Noether 流は全空間について積分しなくても、区間ごとの測定量に意味を持つ。
:::

### 二成分の内部回転が与える保存流

二つの実変位 $u,v$ に同じ慣性と張力を与え、局所的な復元エネルギーを $U(u^2+v^2)$ とする。一様弦上では
$$
{\cal L}=\tfrac{\rho}{2}(u_t^2+v_t^2)-\tfrac{\mathcal T}{2}(u_x^2+v_x^2)-U(u^2+v^2).
$$
媒質上の位置 $x$ を固定して $(u,v)$ の平面だけを回す。$\delta u=-\epsilon v,\delta v=\epsilon u$ だから保存流は
$$
j^0=\rho(uv_t-vu_t),\qquad j^x=-\mathcal T(uv_x-vu_x).
\label{eq-am-field-internal-rotation}
$$
事実 $\partial_tj^0+\partial_xj^x=u(\rho v_{tt}-\mathcal Tv_{xx})-v(\rho u_{tt}-\mathcal Tu_{xx})$。復元力の二項はともに $-2uvU'$ で打ち消される。複素場 $\Phi=(u+iv)/\sqrt2$ と書くと変換は $\Phi\mapsto e^{i\epsilon}\Phi$ である。この :term[内部対称性]{en="internal symmetry" yomi="ないぶたいしょうせい"} は空間を回転する操作と区別される。同じ密度で二方向にたわむ一本の弦なら $j^0$ は媒質の軸まわりの角運動量密度とも解釈できる。$u,v$ が別の内部自由度であれば、同じ式が機械的な角運動量を直接表すとは限らない。

以下では $u=A\cos\zeta,\ v=B\sin\zeta$ と書き、$A,B$ は符号付きの振幅とする。$B$ の符号を変えると、波の進行方向を保ったまま各点での回転の向きだけが変わる。

:::example[直線偏光]
$\rho=2\,\mathrm{kg/m},\mathcal T=50\,\mathrm N$、局所復元力は零とし、$\zeta=2x-10t$、変位を $u=\frac{1}{100}\cos\zeta,\ v=(0)\sin\zeta$（m）とする。波速 $5\,\mathrm{m/s}$ は材料定数からも位相からも同じ値である。二つの変位の積を微分して反対称な組を取れば $\cos^2\zeta+\sin^2\zeta=1$ が残るので

$$
j^0=0\,\mathrm{kg\,m/s},\qquad j^x=0\,\mathrm{N\,m}.
$$

一方、一周期にわたる平均エネルギー密度は $\langle e\rangle=\rho\omega^2(A^2+B^2)/2=\frac{1}{100}\,\mathrm{J/m}$。変位の回転密度も右向きに運ばれ、$j^x=cj^0$ が成り立つ。全エネルギーが零でなくても、内部回転に対する保存量は零になり得る。
:::

:::example[正の向きの円偏光]
$\rho=2\,\mathrm{kg/m},\mathcal T=50\,\mathrm N$、局所復元力は零とし、$\zeta=2x-10t$、変位を $u=\frac{1}{100}\cos\zeta,\ v=(- \frac{1}{100})\sin\zeta$（m）とする。二つの変位の積を微分して反対称な組を取れば $\cos^2\zeta+\sin^2\zeta=1$ が残るので

$$
j^0=\frac{1}{500}\,\mathrm{kg\,m/s},\qquad j^x=\frac{1}{100}\,\mathrm{N\,m}.
$$

一方、一周期にわたる平均エネルギー密度は $\langle e\rangle=\rho\omega^2(A^2+B^2)/2=\frac{1}{50}\,\mathrm{J/m}$。各位置で変位の先端が円を描き、内部の回転密度は時刻によらない。
:::

:::example[逆の向きの楕円偏光]
$\rho=2\,\mathrm{kg/m},\mathcal T=50\,\mathrm N$、局所復元力は零とし、$\zeta=2x-10t$、変位を $u=\frac{1}{50}\cos\zeta,\ v=(\frac{1}{100})\sin\zeta$（m）とする。二つの変位の積を微分して反対称な組を取れば $\cos^2\zeta+\sin^2\zeta=1$ が残るので

$$
j^0=- \frac{1}{250}\,\mathrm{kg\,m/s},\qquad j^x=- \frac{1}{50}\,\mathrm{N\,m}.
$$

一方、一周期にわたる平均エネルギー密度は $\langle e\rangle=\rho\omega^2(A^2+B^2)/2=\frac{1}{20}\,\mathrm{J/m}$。回転の向きを反転すると電荷に相当する保存量は反転するが、エネルギーは振幅の二乗で決まる。
:::

## エネルギー運動量テンソル

### 並進の向きと運動量の符号

時空の一方向 $\nu$ への並進 $x'^\mu=x^\mu+\epsilon\delta^\mu_\nu$、$\Delta^a=0$ を考える。
$\mathcal L$ が $x^\nu$ に陽に依存しなければ、境界項は不要である。
上の定理の符号では $j^\mu=-T^\mu{}_\nu$ となり、
:term[正準エネルギー運動量テンソル]{en="canonical energy-momentum tensor" yomi="せいじゅんえねるぎーうんどうりょうてんそる"} は

$$
T^\mu{}_\nu=\sum_{a=1}^s\Pi_a^\mu\partial_\nu\phi^a-\delta^\mu_\nu\mathcal L .
\label{eq-am-stress-energy}
$$

直接微分すると混合微分項が連鎖律の対応する項を消し、

$$
\sum_\mu\partial_\mu T^\mu{}_\nu
=-\sum_a\mathcal E_a\partial_\nu\phi^a-
\left(\partial_\nu\mathcal L\right)_{\rm explicit}.
$$

従って、外部媒質や外力に陽な位置・時間依存があればそれが収支の源となる。
ここでは $x^0=t$ としたため成分の次元が異なり、この混合成分行列を転置しても一般に同じにならない。

### 弦におけるエネルギーと波の運動量

時間列の保存則は、単位長さのエネルギー $e$ と、右へ通過する仕事率 $S$ により

$$
e=T^0{}_0=\tfrac12(\rho u_t^2+\mathcal T u_x^2),\qquad
S=T^x{}_0=-\mathcal T u_tu_x,\qquad \partial_t e+\partial_x S=0 .
$$

右向き波 $u=A\cos(kx-\omega t)$ では $u_t$ と $u_x$ の符号が反対なので $\langle S\rangle=\mathcal T A^2\omega k/2=c\langle e\rangle>0$。
$A=1.0\,\mathrm{mm},k=10\,\mathrm{m^{-1}},c=100\,\mathrm{m/s},\mathcal T=100\,\mathrm N$ なら $\omega=1000\,\mathrm{s^{-1}}$、
平均仕事率は $0.50\,\mathrm W$ である。
定在波では両方向の輸送が打ち消し、一周期の平均は零になる。

空間並進で「波形を右にずらす」と同じ位置では $\delta u=-\epsilon u_x$ である。
この能動変換の密度 $g=-T^0{}_x=-\rho u_tu_x$ は右向き波で正で、その流束は $-T^x{}_x=e$。
従って $\partial_tg+\partial_x e=0$。
この $g$ は、波形を媒質の座標に対してずらす並進の :term[擬運動量密度]{en="pseudomomentum density" yomi="ぎうんどうりょうみつど"} である。弦の物質点の横運動量密度は $\rho u_t$ であり、こちらは変位そのものを一様にずらす変換から得られる。
$S=c^2g$ はこの一様弦での関係である。

### 保存流の不定性

$U^{\lambda\mu}{}_\nu=-U^{\mu\lambda}{}_\nu$ として $T^\mu{}_\nu$ に $\sum_\lambda\partial_\lambda U^{\lambda\mu}{}_\nu$ を足しても、
二重発散は対称な二階微分と反対称な係数の和なので零になる。
従って局所収支の形は同じである。
全量の差は空間境界項になり、境界で消えるときに限り全量も一致する。
時空計量に関する作用の変分から対称なテンソルを定める方法もある。
ここで求めた正準テンソルは並進に対応する。

### 真空電磁場でゲージに依存しないエネルギーを取り出す

[電磁気学](/physics/electromagnetism/) の SI ポテンシャルを使い、$\vb E=-\partial_t\vb A-\nabla\Phi$、$\vb B=\nabla\times\vb A$ とする。源のない場合の密度は $\mathcal L=\epsilon_0|\vb E|^2/2-|\vb B|^2/(2\mu_0)$ である。$\partial_t\vb A$ に対する共役密度は $-\epsilon_0\vb E$ なので、時間並進に対応する式を途中まで計算すると
$$
T^0{}_0=u_{\rm em}+\epsilon_0\vb E\cdot\nabla\Phi,\qquad
u_{\rm em}=\frac{\epsilon_0}{2}|\vb E|^2+\frac{1}{2\mu_0}|\vb B|^2 .
$$
真空では $\nabla\cdot\vb E=0$ だから余分な項は $\nabla\cdot(\epsilon_0\Phi\vb E)$ に等しい。従って $U^{i0}{}_0=-\epsilon_0\Phi E_i$ と選べば、その発散を打ち消せる。

空間部分も計算できる。$\partial_i\Phi$ に関する微分は $-\epsilon_0 E_i$、$\partial_i A_j$ に関する微分は $-\sum_k\varepsilon_{kij} B_k/\mu_0$ である。ここで $\varepsilon_{kij}$ は右手系の完全反対称記号である。これらを $T^i{}_0$ に入れると、空間流束の列は
$$
\vb S_{\rm em}+\frac{\nabla\Phi\times\vb B}{\mu_0}-\epsilon_0\vb E\dot\Phi,\qquad
\vb S_{\rm em}=\frac{\vb E\times\vb B}{\mu_0}.
$$
真空の Ampère–Maxwell の式 $\nabla\times\vb B/\mu_0=\epsilon_0\partial_t\vb E$ により、余分な流束は $\nabla\times(\Phi\vb B)/\mu_0-\partial_t(\epsilon_0\Phi\vb E)$ となる。後者は $U^{0i}{}_0=-U^{i0}{}_0$ の寄与で消え、前者は $U^{ji}{}_0=-\sum_k\varepsilon_{ijk}\Phi B_k/\mu_0$ の空間発散で消える。残った密度と流束はどちらも $\Phi,\vb A$ のゲージを変えても同じである。

例えば真空の平面波では $\vb B=\hat{\vb k}\times\vb E/c$、$c^{-2}=\epsilon_0\mu_0$ だから、瞬間の密度は $u_{\rm em}=\epsilon_0|\vb E|^2$、流束は $c u_{\rm em}\hat{\vb k}$ となる。保存流の改善が、実際に電場と磁場で測れるエネルギー輸送を与える例である。

### 波の運動量収支を境界で読む

固定端では $u_t=0$ なのでエネルギー流は零である。
それでも傾きが零とは限らず、空間並進に対応する流束 $-T^x{}_x=e$ は正になりうる。
支持具は動かないので仕事をしないが、波形の並進対称性を破る。
エネルギー保存から空間並進の保存まで推測してはならない。
自由端では $u_x=0$ で同じくエネルギー流は零だが、こちらでは端が上下に動ける。

非一様弦で $\rho(x)$ と $\mathcal T(x)$ が時間に依存しないとき、時間列は保存する。
空間列には源が残る。
$g=-\rho u_tu_x$ として式を展開すると

$$
\partial_t g+\partial_x e
=\tfrac12\rho'(x)u_t^2-\tfrac12\mathcal T'(x)u_x^2 .
$$

これは媒質の形によって並進の対称性が失われる程度を定量化する。
$\mathcal T$ だけ滑らかに変化する場合でも、エネルギーを保ったまま波の運動量は変わりうる。

### 回転に対する保存流

三次元での単一実スカラー場を考え、

$$
\mathcal L=\frac{\rho}{2}\phi_t^2-\frac{\kappa}{2}|\nabla\phi|^2-U(\phi)
$$

とする。
$\rho,\kappa$ は正の等方的定数である。
$z$ 軸まわりの座標回転なら $(\xi^x,\xi^y,\xi^z)=(-y,x,0)$、$\xi^0=0,\Delta=0$ である。
場はスカラーなので対応する点で成分自体の回転はない。
固定位置での変化は $\eta=y\phi_x-x\phi_y$、従って

$$
j^0=\rho\phi_t(y\phi_x-x\phi_y)
=xg_y-yg_x,\qquad g_i=-\rho\phi_t\partial_i\phi .
$$

質点で位置と運動量の外積を取ったのと同じ構造である。
ベクトル場の場合には成分そのものも回転し、$\Delta^a$ の項が流れに加わる。
この違いは空間の添字と場の成分番号を区別してはじめて見える。

### 分散する弦ではどの速度でエネルギーが運ばれるか

床ばねの密度を $\kappa u^2/2$ とした一様弦で、全エネルギー密度は

$$
e=\tfrac12(\rho u_t^2+\mathcal T u_x^2+\kappa u^2).
$$

流束は空間微分に対する共役密度で決まるため $S=-\mathcal T u_tu_x$ のままである。
$u=A\cos(kx-\omega t)$ で $\rho\omega^2=\mathcal T k^2+\kappa$ を使うと、

$$
\begin{aligned}
\langle e\rangle&=\tfrac14 A^2(\rho\omega^2+\mathcal T k^2+\kappa)
=\tfrac12\rho A^2\omega^2,\\
\langle S\rangle&=\tfrac12\mathcal T A^2 k\omega .
\end{aligned}
$$

従ってエネルギー流とエネルギー密度の比は

$$
\frac{\langle S\rangle}{\langle e\rangle}
=\frac{\mathcal T k}{\rho\omega}
=\frac{\dd\omega}{\dd k}.
$$

これは群速度であり、$\kappa>0$ のとき位相速度 $\omega/k$ と異なる。
$k=0$ では全点が同位相で上下し、有限なエネルギーを持っていても流束は零である。

:::example[離散化の間隔を決める]
$\rho=1\,\mathrm{kg/m},\mathcal T=4\,\mathrm N,\kappa=9\,\mathrm{N/m^2}$、$k=2\,\mathrm{m^{-1}}$ なら $\omega=5\,\mathrm{s^{-1}}$。
位相速度は $2.5\,\mathrm{m/s}$、群速度は $1.6\,\mathrm{m/s}$ である。
連続弦を格子間隔 $a=0.05\,\mathrm m$ で近似すると $ka=0.1$。
差分の弾性項の正確な比は $[2\sin(ka/2)/(ka)]^2$ で、連続値からの相対差は約 $(ka)^2/12=8.3\times10^{-4}$。
物理的な床ばねによる分散と、格子を粗くしたことによる数値的な分散を区別できる。
:::

### 弾性体で応力が運ぶエネルギー

基準配置の小領域に注目すると、表面を横切って何が運ばれるかは応力と境界速度から決まる。運動方程式を $\dot u_i$ 倍して $i$ について和を取り、積の微分を使えば
$$
\partial_t e+\sum_j\partial_j S_j=\sum_i f_i\dot u_i,\quad
e=\frac\rho2\sum_i\dot u_i^2+W,\qquad S_j=-\sum_i\sigma_{ij}\dot u_i .
$$
実際、$\dot W=\sum_{i,j}\sigma_{ij}\partial_j\dot u_i$ なので内部の応力の仕事が蓄積エネルギーの時間微分になる。面の外向きエネルギー流は $\vb n\cdot\vb S=-\vb t\cdot\dot{\vb u}$ である。外の物体が内側に正の仕事をすると、外向き流束は負である。$\sigma$ は媒質の力を数えるのに対し、$T^\mu{}_\nu$ の空間並進は基準配置上で変位の模様をずらす対称性に対応する。例えば均一な静的引張りでは応力が有限でもエネルギーの流れは零である。

縦または横の平面進行波で変位振幅を $A$、速度を $c$ と書けば、ある伝播方向に沿う有効剛性は $C=\rho c^2$。傾きの振幅が $kA$ なので弾性エネルギーの時間平均は $C k^2 A^2/4$、運動部分は $\rho\omega^2A^2/4$。分散関係 $\omega=ck$ で両者が等しい。従って平均強度は $\langle S\rangle=\rho c\omega^2A^2/2$。例えば $\rho=1000\,\mathrm{kg/m^3},c=1000\,\mathrm{m/s},\omega=1000\,\mathrm{s^{-1}},A=10^{-6}\,\mathrm m$ なら $0.50\,\mathrm{W/m^2}$ になる。二つの偏光が重なるときはまず変位と応力を足し、その積を平均する。干渉項が消える条件を確かめずに各波の強度を加えることはできない。

## 場の Hamilton 形式

### どの微分を Legendre 変換するか

時間発展の初期値を指定するには、一つの時刻で場とその速度を指定する。
そのため $\phi^a_t$ だけについて :ref[def-am-hamiltonian] の変換を行う。
空間微分と陽に現れる位置 $\vb x$ はそのまま引数に残し、

$$
\pi_a=\Pi_a^0,\qquad
\mathcal H(\phi,\pi,\nabla\phi,\vb x,t)=\sum_a\pi_a\phi^a_t-\mathcal L,\qquad
H=\int_D\mathcal H\,\dd^dx
$$

と定義する。
速度 Hessian が正則で $\phi_t$ が局所的に消去できると仮定した。
$\mathcal H=T^0{}_0$ である。
弦では $\pi=\rho u_t$、$\mathcal H=\pi^2/(2\rho)+\mathcal T u_x^2/2$。

### 汎関数微分と等時刻括弧

$F[\phi,\pi]$ の :term[汎関数微分]{en="functional derivative" yomi="はんかんすうびぶん"} は、時刻を固定した第一変分を

$$
\delta F=\int_D\sum_a\left(\frac{\delta F}{\delta\phi^a(\vb x)}\delta\phi^a(\vb x)+
\frac{\delta F}{\delta\pi_a(\vb x)}\delta\pi_a(\vb x)\right)\dd^dx
$$

と書く係数である。
境界項を消す許容変分を取る。
$H$ の微分では $\delta(\partial_i\phi^a)$ の項を部分積分して

$$
\frac{\delta H}{\delta\phi^a}=\mathcal H_{\phi^a}-\sum_i\partial_i\mathcal H_{\partial_i\phi^a},\qquad
\frac{\delta H}{\delta\pi_a}=\mathcal H_{\pi_a}.
$$

$\int[\int_D\sum_a\pi_a\dot\phi^a\,\dd^dx-H]\dd t$ を独立な $\phi,\pi$ で変分し、時間端の $\delta\phi$ を零とすれば

$$
\dot\phi^a=\frac{\delta H}{\delta\pi_a},\qquad
\dot\pi_a=-\frac{\delta H}{\delta\phi^a}.
\label{eq-am-field-hamilton}
$$

Legendre 変換により $\mathcal H_{\phi^a}=-\mathcal L_{\phi^a}$、
$\mathcal H_{\partial_i\phi^a}=-\Pi_a^i$ だから第二式は場の Euler–Lagrange 方程式に一致する。
弦では $\dot u=\pi/\rho,\ \dot\pi=\mathcal T u_{xx}$ が再現される。

固定時刻での :ref[def-am-poisson] を場に移すと

$$
\{F,G\}=\int_D\sum_a\left(
\frac{\delta F}{\delta\phi^a}\frac{\delta G}{\delta\pi_a}-
\frac{\delta F}{\delta\pi_a}\frac{\delta G}{\delta\phi^a}
\right)\dd^dx .
$$

試験関数 $f,g$ と場の積分に適用すると $\{\int f\phi^a,\int g\pi_b\}=\delta^a_b\int fg$。
この性質を分布で $\{\phi^a(\vb x),\pi_b(\vb y)\}=\delta^a_b\delta^{(d)}(\vb x-\vb y)$ と書く。
同種同士の括弧は零である。
例えば周期境界または全空間で境界項が消える場合、$P=-\int\pi u_x\dd x$ なら $\{u(x),P\}=-u_x(x)$ であり、空間並進の符号を確かめられる。

### 離散系の括弧がデルタ関数になるまで

離散化した弦の格子点で共役運動量は $p_j=\rho a\dot u_j$ である。
一点の質量に対しては $\{u_j,p_k\}=\delta_{jk}$。
一方、運動量**密度**の標本は $\pi_j=p_j/a$ だから

$$
\{u_j,\pi_k\}=\frac{\delta_{jk}}a,\qquad
\sum_k a\,\frac{\delta_{jk}}a f_k=f_j .
$$

最後の和が連続極限で $\int\delta(x-y)f(y)\dd y=f(x)$ になる。
$1/a$ を含める理由は格子の自由度が増えることだけでなく、点ごとの運動量と長さあたりの密度が異なるためである。
弦で $\delta(x-y)$ の単位は $\mathrm{m^{-1}}$ となる。

### 基準モードの正準規格化

両端固定弦で $f_n(x)=\sqrt{2/l}\sin(n\pi x/l)$ と置く。
三角関数の積を和に直して積分すれば $\int_0^l f_n f_m\dd x=\delta_{nm}$。
$Q_n=\int f_n u\dd x,\ P_n=\int f_n\pi\dd x$ とすると等時刻括弧から

$$
\{Q_n,P_m\}=\int_0^l f_n(x)f_m(x)\dd x=\delta_{nm} .
$$

滑らかな固定端変位を Fourier 正弦級数 $u=\sum_n Q_nf_n$ に展開する。
運動量についても同じ展開を用い、部分積分と $f_n''=-(n\pi/l)^2 f_n$ を使うと

$$
H=\sum_{n=1}^\infty\left[\frac{P_n^2}{2\rho}
+\frac{\mathcal T}{2}\left(\frac{n\pi}l\right)^2 Q_n^2\right].
$$

有限エネルギーの場では各項の和が収束する。
例えば $Q_2,P_2$ のみを初期に非零にすると、その後も第二モードだけの調和振動が続く。
弦の全自由度は無限個でも、線形性によって初期に励起されたモードを一つずつ追えばよい。

### 速度が逆に解けない場の例

電磁場をスカラーポテンシャル $\Phi$ とベクトルポテンシャル $\vb A$ で書くと、真空中の SI 密度は

$$
\mathcal L=\frac{\epsilon_0}{2}|-\nabla\Phi-\partial_t\vb A|^2
-\frac{1}{2\mu_0}|\nabla\times\vb A|^2 .
$$

$\partial_t\Phi$ が現れないので $\pi_\Phi=0$。
一方 $\pi_{A_i}=-\epsilon_0 E_i$ である。
$\Phi$ の Euler–Lagrange 方程式は $\sum_i\partial_i(\epsilon_0 E_i)=0$、すなわち電荷のない Gauss の法則となる。
従って四成分を独立な四つの通常の正準対とみなすことはできない。
このように速度の逆変換ができないときは、拘束を時間発展と両立させる Hamilton 形式が必要になる。
その一般的な構成は :cite[onuki-am] に譲る。
[電磁気学](/physics/electromagnetism/)

### Hamiltonian を微分できる境界データ

生成子を場の括弧に入れるには、領域内の汎関数微分が存在するだけでなく、許容変分に対して境界項も消えなければならない。
例えば区間 $[0,l]$ の弦で

$$
\delta H=\int_0^l
\left(\frac{\pi}{\rho}\delta\pi-\mathcal T u_{xx}\delta u\right)\dd x
+[\mathcal T u_x\delta u]_0^l .
$$

指定した固定端なら $\delta u=0$ で消える。
端の変位を許し、端に既知の一定力 $f$ を加える場合には $H_{\rm total}=H-fu(l)$ と置く。
左端は固定すると、右端の変分は $(\mathcal T u_x(l)-f)\delta u(l)$。
つり合い条件 $\mathcal T u_x(l)=f$ がちょうど生成子を微分可能にする。
外から端の位置を時間依存に指定する場合には、許容な場の時間発展そのものが境界を動かすため、エネルギーには支持具の仕事も数える。

### 局所密度から生成子を組み立てる

周期境界の二成分弦の運動量密度を $\pi_u,\pi_v$ とする。
内部回転の全保存量は

$$
J=\int_0^l(u\pi_v-v\pi_u)\dd x .
$$

汎関数微分から

$$
\begin{aligned}
\{u(x),J\}&=-v(x),& \{v(x),J\}&=u(x),\\
\{\pi_u(x),J\}&=-\pi_v(x),& \{\pi_v(x),J\}&=\pi_u(x).
\end{aligned}
$$

運動量まで同じように回転することが分かる。
同じ $\rho,\mathcal T$ を両成分に与えた Hamiltonian では $\pi_u^2+\pi_v^2$ と $u_x^2+v_x^2$ が不変なので $\{H,J\}=0$。
局所流から出発して積分した量が、場の Hamilton 形式で実際に対称性を生成することを確認できた。
二成分で異なる床ばねを付ければこの内部回転は破れ、同じ $J$ の時間微分には二つのばね力の差が残る。

### 電磁場の正準運動量と Gauss の拘束

波動方程式には二階の時間微分が現れるが、場の成分すべてに独立な運動量があるとは限らない。SI 単位で $\vb E=-\partial_t\vb A-\nabla\Phi,\ \vb B=\nabla\times\vb A$ と置き、電荷・電流密度を $\rho_e,\vb j$ とする。外部源を指定した作用密度は
$$
{\cal L}=\frac{\epsilon_0}{2}|\vb E|^2-\frac{1}{2\mu_0}|\vb B|^2-\rho_e\Phi+\vb j\cdot\vb A .
$$
場の共役運動量は $\boldsymbol\pi=-\epsilon_0\vb E,\ \pi_\Phi=0$ である。$\dot{\vb A}=\boldsymbol\pi/\epsilon_0-\nabla\Phi$ は得られるが $\dot\Phi$ は得られない。境界変分と表面項が消えるよう周期境界または十分速い減衰を仮定し、$\int-\boldsymbol\pi\cdot\nabla\Phi$ を部分積分すれば
$$
H_T=\int\left[\frac{|\boldsymbol\pi|^2}{2\epsilon_0}+
\frac{|\nabla\times\vb A|^2}{2\mu_0}+\Phi(\nabla\cdot\boldsymbol\pi+\rho_e)-
\vb j\cdot\vb A+w\pi_\Phi\right]\dd^3x .
\label{eq-am-maxwell-hamiltonian}
$$
$w$ は未定の速度である。$\dot\pi_\Phi=-\delta H_T/\delta\Phi=0$ より $\nabla\cdot\boldsymbol\pi=-\rho_e$、すなわち Gauss の式が出る。これは初期値に課す条件である。周期領域では発散の積分が零になるため、指定できる全電荷も零に限られる。$\vb A$ については $\int\vb B\cdot(\nabla\times\delta\vb A)=\int\delta\vb A\cdot(\nabla\times\vb B)$ より $\dot{\boldsymbol\pi}=\vb j-\nabla\times\vb B/\mu_0$ を得る。この式の発散は $\nabla\cdot\dot{\boldsymbol\pi}=\nabla\cdot\vb j$ だから、拘束が保存するための条件は $\dot\rho_e+\nabla\cdot\vb j=0$。外部源を独立に好きな関数として指定してよいわけではない。

ゲージ関数 $\chi$ に対して $\delta\vb A=\nabla\chi,\delta\Phi=-\dot\chi$ を起こす等時刻の生成子は $G=\int[-\chi(\nabla\cdot\boldsymbol\pi+\rho_e)-\dot\chi\,\pi_\Phi]\dd^3x$。$\chi$ は境界で消えるとする。部分積分から $\delta G/\delta\boldsymbol\pi=\nabla\chi,\ \delta G/\delta\pi_\Phi=-\dot\chi$ を確かめられる。$\vb E,\vb B$ は変わらず、源との結合の変化は電荷保存のもとで境界項になる。真空の非零 Fourier 波数では $\vb A$ の縦成分をゲージで消し、$\boldsymbol\pi$ の縦成分は Gauss 条件で消える。$\Phi,\pi_\Phi$ の組にも独立な波はない。従って自由な波の偏光は横の二つである。一定場など波数零のモードと境界に残る電束は、形状と境界条件を別に指定して扱う。

## 演習問題

::::exercise{level=basic}
$\mathcal T=40\,\mathrm N,\rho=0.10\,\mathrm{kg/m},l=2\,\mathrm m$ の両端固定弦の基本周波数を求めよ。

:::solution
$c=\sqrt{40/0.10}=20\,\mathrm{m/s}$、基本波長は $2l=4\,\mathrm m$ なので $f_1=5\,\mathrm{Hz}$。
:::
::::

::::exercise{level=basic}
左端固定、右端自由な長さ $l$ の弦の許容波数を求めよ。

:::solution
左端の条件から空間解は $\sin kx$。
右端で $k\cos kl=0$ より $k=(n+1/2)\pi/l,\ n=0,1,\ldots$。
$k=0$ では線形な静的解に両端条件を課すと零に限られる。
:::
::::

::::exercise{level=basic}
$\rho,\mathcal T$ 一定の弦に長さあたりの外力 $f(x,t)$ がある。
密度にどの項を加えるか。
またエネルギー供給率を求めよ。

:::solution
$+fu$ を加えると $\rho u_{tt}=\mathcal T u_{xx}+f$。
両辺に $u_t$ を掛けて積の微分に直せば、弦自身のエネルギーについて $\partial_te+\partial_xS=fu_t$ である。
:::
::::

::::exercise{level=basic}
$u=A\sin kx\cos\omega t$ のエネルギー流の一周期平均を求めよ。

:::solution
$S=\mathcal T A^2\omega k\sin kx\cos kx\sin\omega t\cos\omega t$。
時間因子は $\sin(2\omega t)/2$ なので平均は零である。
:::
::::

::::exercise{level=advanced}
$\mathcal T(x)>0$ が変化する弦の方程式と、時間不変の張力でのエネルギー収支を導け。

:::solution
$\Pi^x=-\mathcal T(x)u_x$ から $\rho u_{tt}=\partial_x(\mathcal T u_x)$。
$e=\rho u_t^2/2+\mathcal T u_x^2/2$ を微分し、方程式に $u_t$ を掛けたものを代入すると $e_t=\partial_x(\mathcal T u_xu_t)$。
$\mathcal T'$ は空間微分の中に含まれる。
:::
::::

::::exercise{level=advanced}
右端に点質量 $M$ を付け、端の横運動を自由にする。
境界条件から全エネルギーが保存することを示せ。
左端は固定する。

:::solution
端の作用に $M u_t(l)^2/2$ を加える。
時間について部分積分すると境界変分は $[-M u_{tt}(l)-\mathcal T u_x(l)]\delta u(l)$ だから $M u_{tt}(l)=-\mathcal T u_x(l)$。
端の運動エネルギーの微分は $-\mathcal T u_x(l)u_t(l)=S(l)$。
弦のエネルギー変化は $S(0)-S(l)=-S(l)$ であり和は一定である。
:::
::::

::::exercise{level=advanced}
平方可積分な周期場で $P=-\int\pi u_x\dd x$ が $u,\pi$ の両方を同じ向きにずらすことを示せ。
場は滑らかとする。

:::solution
変分して $\delta P=\int[-u_x\delta\pi+\pi_x\delta u]\dd x$。
周期条件で部分積分の端点は相殺される。
従って $\{u,P\}=-u_x$、$\{\pi,P\}=-\pi_x$ であり、両方とも $x$ で評価した値が $x-\epsilon$ の値に移る。
:::
::::

::::exercise{level=advanced}
固定端弦で $u=Q(t)\sin(n\pi x/l)$ の一モードだけを残した Hamiltonian を求めよ。

:::solution
$\sin^2,\cos^2$ の積分はいずれも $l/2$。
よって $L=(\rho l/4)\dot Q^2-(\mathcal Tl/4)(n\pi/l)^2Q^2$。
有効質量 $m_n=\rho l/2$、運動量 $P=m_n\dot Q$ とすれば $H=P^2/(2m_n)+m_n\omega_n^2Q^2/2$、
$\omega_n=(n\pi/l)\sqrt{\mathcal T/\rho}$ である。
:::
::::

::::exercise{level=advanced}
時間一定の密度 $\rho$ と、一様だが時間で変化する張力 $\mathcal T(t)$ をもつ弦で、エネルギー密度と空間並進の密度の収支を求めよ。
周期境界を仮定する。

:::solution
局所方程式は $\rho u_{tt}=\mathcal T(t)u_{xx}$。
$e=\rho u_t^2/2+\mathcal T(t)u_x^2/2$ の微分には陽な時間変化が残るため、

$$
\partial_t e+\partial_x(-\mathcal T u_tu_x)=\dot{\mathcal T}\,u_x^2/2 .
$$

一方 $g=-\rho u_tu_x$ については $\partial_t g+\partial_x e=0$ が成り立つ。
周期境界で積分すれば、全エネルギーには外部装置からの仕事が加わるが、波形の空間並進を生成する全運動量は一定である。
:::
::::

::::exercise{level=advanced}
相互作用する二つの変位場 $u,v$ を、同じ自由弦密度の和に $-\kappa(u-v)^2/2$ を加えて記述する。
保存される内部並進と、二つの基準場の分散関係を求めよ。

:::solution
$u$ と $v$ を同じ定数だけ増やせば密度は不変なので、保存密度は $\rho(u_t+v_t)$、流束は $-\mathcal T(u_x+v_x)$ である。
一つずつの場の方程式は

$$
\rho u_{tt}=\mathcal T u_{xx}-\kappa(u-v),\qquad
\rho v_{tt}=\mathcal T v_{xx}+\kappa(u-v).
$$

$w_+=(u+v)/\sqrt2,\ w_-=(u-v)/\sqrt2$ に替え、$e^{i(kx-\omega t)}$ を入れると $\omega_+^2=\mathcal T k^2/\rho,\ \omega_-^2=(\mathcal T k^2+2\kappa)/\rho$。
空間一様な同位相変位は自由だが、逆位相変位では局所ばねが伸びる。
:::
::::

::::exercise{level=advanced}
周期領域の真空電磁場について、上の Hamilton 方程式からエネルギー密度と局所エネルギー流を導け。外部電流を戻したときの場の仕事率の符号も決めよ。

:::solution
$\epsilon_0\dot{\vb E}=(\nabla\times\vb B)/\mu_0-\vb j$ を $\vb E$ と内積し、$\dot{\vb B}=-\nabla\times\vb E$ を $\vb B/\mu_0$ と内積する。和を取ると $\partial_t(\epsilon_0 E^2/2+B^2/(2\mu_0))=[\vb E\cdot(\nabla\times\vb B)-\vb B\cdot(\nabla\times\vb E)]/\mu_0-\vb j\cdot\vb E$。角括弧は $-\nabla\cdot(\vb E\times\vb B)$ なので流束は $\vb E\times\vb B/\mu_0$。周期境界では領域全体の流出が相殺し、物質が単位体積・単位時間あたり $\vb j\cdot\vb E$ を受け取るだけ場のエネルギーが減る。
:::
::::

## 確認テスト

::::quiz{answer="30" tol="1e-6"}
$\mathcal T=90\,\mathrm N,\rho=0.10\,\mathrm{kg/m}$ で波の速さを m/s で答えよ。

:::explanation
$\sqrt{\mathcal T/\rho}=\sqrt{900}=30\,\mathrm{m/s}$ である。
:::
::::

::::quiz
外力も端質量もない自由端の自然境界条件はどれか。
張力は正とする。

- [ ] $u=0$
- [x] $u_x=0$
- [ ] $u_t=0$

:::explanation
端の変位の変分が任意であるため、その係数 $-\mathcal T u_x$ が零になる。
:::
::::

::::quiz{answer="4" tol="1e-6"}
$\rho=2\,\mathrm{kg/m},u_t=2\,\mathrm{m/s}$ で弦の共役運動量密度の値を SI 単位で答えよ。

:::explanation
$\pi=\partial\mathcal L/\partial u_t=\rho u_t=4\,\mathrm{kg/s}$ である。
:::
::::

::::quiz
局所的な保存則から有限領域の全量が一定となる条件はどれか。

- [ ] 場が時間に依存する
- [ ] 密度が全点で正である
- [x] 境界を通る正味の流束が零である

:::explanation
発散定理により全量の時間微分は外向き流束の面積分の負に等しいので、その積分が消えることが条件である。
:::
::::

::::quiz
二本の弦を点質量のない結び目でつなぎ、結び目に横方向の外力が働かないとする。
横変位が連続であるほかに必要な条件はどれか。

- [ ] 密度が左右で等しい
- [x] $\mathcal T u_x$ が左右で連続する
- [ ] 波速が左右で等しい

:::explanation
左右の端点項で結び目の変分は共通である。
外力も集中質量もなければ、その係数である横張力の差が零になる。
:::
::::

::::quiz{answer="0.75" tol="1e-6"}
波動インピーダンスが右で左の三倍となる界面に波が入射する。
透過エネルギー流の入射エネルギー流に対する割合を答えよ。

:::explanation
透過振幅比は $2Z_1/(Z_1+3Z_1)=1/2$。
エネルギー流は $Z|A|^2$ に比例するため、透過の割合は $3(1/2)^2=0.75$ である。
:::
::::
